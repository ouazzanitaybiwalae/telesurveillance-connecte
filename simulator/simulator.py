#!/usr/bin/env python3
"""
Simulateur de dispositifs connectes - CTSMC (Centre de TeleSuivi des Maladies Chroniques)
Projet 8 : systeme de telesurveillance connecte, suivi du diabete (type 1 et type 2).

Simule la PASSERELLE / SMARTPHONE du patient (client MQTT) qui relaie les mesures de :
    - glucometre      IEEE 11073-10417   -> ctmc/patients/{pseudo_id}/glycemie
    - tensiometre     IEEE 11073-10407   -> ctmc/patients/{pseudo_id}/tension
    - thermometre     IEEE 11073-10408   -> ctmc/patients/{pseudo_id}/temperature
    - statut (LWT)                       -> ctmc/devices/{gateway_id}/status  (retained)

METHODE DE GENERATION DES DONNEES (a recopier dans le rapport)
--------------------------------------------------------------
Toutes les donnees sont FICTIVES. Identifiants pseudonymes : P-7f3a9c, P-2b81d4, P-c04e57.
Pour chaque grandeur, a chaque pas de temps :
    x[t+1] = x[t] + theta * (baseline - x[t]) + sigma * N(0,1)     marche aleatoire avec retour
                                                                    a la moyenne (Ornstein-Uhlenbeck
                                                                    discretise)
    mesure = x[t+1] (+ pic post-prandial pour la glycemie) + bruit * N(0,1)
- Glycemie : un pic post-prandial (gaussienne, amplitude 35-45 mg/dL) revient periodiquement.
- Anomalie : la baseline glisse vers une valeur cible (rampe -> plateau -> rampe de retour).
  Hypoglycemie : la glycemie chute et le pouls monte (tachycardie reactionnelle).
- Temps accelere : un pas de simulation (--interval, 5 s par defaut) represente une prise de
  mesure ; les durees de la simulation ne sont pas des durees cliniques reelles.
- Degradations : panne brutale de la passerelle (LWT), valeurs implausibles, retrait du
  consentement, perte de paquets, broker coupe (tampon local + renvoi).

UNITES : le glucometre 11073-10417 publie en mg/dL (MDC_DIM_MILLI_G_PER_DL).
Les seuils CTSMC sont en g/L : g/L = mg/dL / 100  (0,70 g/L = 70 mg/dL).

CODES MDC : `mdc_code` est le code 32 bits decimal (partition*65536 + code terme), ex. glycemie
MDC_CONC_GLU_CAPILLARY_WHOLEBLOOD = 160184, TA systolique 150021, diastolique 150022,
pouls non invasif 149546, temperature corporelle 150364. Voir HL7 PHD IG + RTMMS (NIST).

FICHIER DE VERITE TERRAIN : results/events.csv (debut/fin d'anomalie, messages perdus,
valeurs implausibles injectees...). Sert a calculer taux de detection, faux positifs, latence.

USAGE
    set MQTT_PASSWORD=...                (cmd)    /   $env:MQTT_PASSWORD="..."  (PowerShell)
    python simulator.py --scenario normal
    python simulator.py --scenario hypoglycemia
    python simulator.py --scenario hyperglycemia
    python simulator.py --scenario hyperglycemia_critical
    python simulator.py --scenario sensor_offline
    python simulator.py --scenario implausible
    python simulator.py --scenario no_consent
    python simulator.py --scenario packet_loss
    python simulator.py --scenario broker_down      (puis docker stop mosquitto / docker start mosquitto)
    python simulator.py --scenario hypoglycemia --dry-run      (aucun broker requis)
    TLS : ajouter  --tls --cafile mosquitto/certs/ca.crt --port 8883   (etape securite)
"""

import argparse
import csv
import getpass
import json
import math
import os
import random
import sys
import time
from collections import deque
from datetime import datetime, timezone
from pathlib import Path

# ----------------------------------------------------------------------------
# Modele de donnees
# ----------------------------------------------------------------------------
ROOT = "ctmc"

# Grandeurs mesurees. mdc_code = code MDC 32 bits en decimal (partition x 65536 + code terme),
# comme dans le mapping FHIR PHD (system urn:iso:std:iso:11073:10101). Verifies sur les exemples
# de l'IG HL7 PHD. Les unites MDC (unit) gardent leur nom ; ucum = equivalent UCUM.
OBS = {
    "glucose":     dict(mdc_name="MDC_CONC_GLU_CAPILLARY_WHOLEBLOOD", mdc_code=160184,
                        unit="MDC_DIM_MILLI_G_PER_DL", ucum="mg/dL", theta=0.20, sigma=1.5, noise=1.0,
                        lo=20, hi=600, dec=0),
    "bp_sys":      dict(mdc_name="MDC_PRESS_BLD_NONINV_SYS", mdc_code=150021,
                        unit="MDC_DIM_MMHG", ucum="mm[Hg]", theta=0.10, sigma=1.5, noise=1.0,
                        lo=60, hi=250, dec=0),
    "bp_dia":      dict(mdc_name="MDC_PRESS_BLD_NONINV_DIA", mdc_code=150022,
                        unit="MDC_DIM_MMHG", ucum="mm[Hg]", theta=0.10, sigma=1.0, noise=0.8,
                        lo=30, hi=150, dec=0),
    "pulse":       dict(mdc_name="MDC_PULS_RATE_NON_INV", mdc_code=149546,
                        unit="MDC_DIM_BEAT_PER_MIN", ucum="/min", theta=0.15, sigma=1.2, noise=0.8,
                        lo=30, hi=220, dec=0),
    "temperature": dict(mdc_name="MDC_TEMP_BODY", mdc_code=150364,
                        unit="MDC_DIM_DEGC", ucum="Cel", theta=0.05, sigma=0.03, noise=0.03,
                        lo=30, hi=43, dec=1),
}

# Un type de message = un topic = un dispositif. `every` : un message tous les N pas.
GROUPS = {
    "glycemie":    dict(prefix="GLU", spec="IEEE 11073-10417 (glucose meter)",
                        metrics=["glucose"], every=1),
    "tension":     dict(prefix="BPM", spec="IEEE 11073-10407 (blood pressure monitor)",
                        metrics=["bp_sys", "bp_dia", "pulse"], every=6),
    "temperature": dict(prefix="THM", spec="IEEE 11073-10408 (thermometer)",
                        metrics=["temperature"], every=6),
}

PATIENTS = {
    "P-7f3a9c": dict(profile="T1", idx=1, meal_amp=45,
                     base=dict(glucose=105, bp_sys=118, bp_dia=76, pulse=72, temperature=36.8)),
    "P-2b81d4": dict(profile="T2", idx=2, meal_amp=40,
                     base=dict(glucose=130, bp_sys=135, bp_dia=85, pulse=78, temperature=36.9)),
    "P-c04e57": dict(profile="T2", idx=3, meal_amp=35,
                     base=dict(glucose=112, bp_sys=128, bp_dia=82, pulse=75, temperature=36.7)),
}

# Scenarios.
#  shift   : la baseline de grandeurs glisse vers des cibles
#  crash   : coupure brutale de la passerelle (declenche le LWT du broker)
#  glitch  : valeurs de glycemie implausibles injectees
#  consent : le patient retire son consentement
#  loss    : perte aleatoire de messages (tous les patients)
#  manual  : rien d'injecte ; on coupe le broker a la main
SCENARIOS = {
    "normal":                 None,
    "hypoglycemia":           dict(kind="shift", patient="P-7f3a9c",
                                   targets=dict(glucose=50, pulse=108)),
    "hyperglycemia":          dict(kind="shift", patient="P-2b81d4", targets=dict(glucose=290)),
    "hyperglycemia_critical": dict(kind="shift", patient="P-2b81d4", targets=dict(glucose=540)),
    "sensor_offline":         dict(kind="crash", patient="P-c04e57"),
    "implausible":            dict(kind="glitch", patient="P-2b81d4"),
    "no_consent":             dict(kind="consent", patient="P-c04e57"),
    "packet_loss":            dict(kind="loss", patient=None),
    "broker_down":            dict(kind="manual", patient=None),
}

MEAL_PERIOD = 150.0   # s, periode des repas simules
MEAL_PEAK = 75.0      # s, instant du pic dans la periode
MEAL_WIDTH = 20.0     # s, largeur du pic
MAX_BUFFER = 1000     # messages gardes en memoire si le broker est injoignable


# ----------------------------------------------------------------------------
# Utilitaires
# ----------------------------------------------------------------------------
def now_iso():
    return datetime.now(timezone.utc).isoformat(timespec="milliseconds").replace("+00:00", "Z")


class EventLog:
    """Journal de verite terrain (CSV)."""

    def __init__(self, path, run_id):
        self.run_id = run_id
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        new = not path.exists()
        self.f = open(path, "a", newline="", encoding="utf-8")
        self.w = csv.writer(self.f)
        if new:
            self.w.writerow(["run_id", "ts", "patient_id", "event", "detail"])

    def log(self, patient_id, event, detail=""):
        self.w.writerow([self.run_id, now_iso(), patient_id, event, detail])
        self.f.flush()

    def close(self):
        self.f.close()


class PatientModel:
    def __init__(self, pid, cfg, rng):
        self.pid, self.cfg, self.rng = pid, cfg, rng
        self.state = dict(cfg["base"])

    def step(self, targets, factor, t):
        """Avance d'un pas ; renvoie les valeurs mesurees (arrondies, bornees)."""
        out = {}
        for name, p in OBS.items():
            baseline = self.cfg["base"][name]
            if name in targets:
                baseline += (targets[name] - baseline) * factor
            x = self.state[name]
            x += p["theta"] * (baseline - x) + p["sigma"] * self.rng.gauss(0, 1)
            self.state[name] = x
            v = x
            if name == "glucose":
                phase = t % MEAL_PERIOD
                bump = self.cfg["meal_amp"] * math.exp(-(((phase - MEAL_PEAK) / MEAL_WIDTH) ** 2))
                v += bump * (1.0 - factor if "glucose" in targets else 1.0)
            v += p["noise"] * self.rng.gauss(0, 1)
            v = max(p["lo"], min(p["hi"], v))
            out[name] = int(round(v)) if p["dec"] == 0 else round(v, p["dec"])
        return out


def ramp_factor(t, start, duration, ramp):
    """0 -> 1 (rampe) -> plateau -> 1 -> 0 (rampe de retour)."""
    end = start + duration
    if t < start:
        return 0.0
    if t < start + ramp:
        return (t - start) / ramp
    if t < end:
        return 1.0
    if t < end + ramp:
        return 1.0 - (t - end) / ramp
    return 0.0


def build_message(pid, group, values, consent=True):
    g = GROUPS[group]
    n = PATIENTS[pid]["idx"]
    observations = []
    for m in g["metrics"]:
        o = OBS[m]
        observations.append({
            "mdc_name": o["mdc_name"],
            "mdc_code": o["mdc_code"],
            "value": values[m],
            "unit": o["unit"],
            "unit_ucum": o["ucum"],
            "quality": "valid",
        })
    return {
        "device_id": f"{g['prefix']}-{n:03d}",
        "device_spec": g["spec"],
        "patient_pseudo_id": pid,
        "profile": PATIENTS[pid]["profile"],
        "consent": consent,
        "timestamp": now_iso(),
        "observations": observations,
    }


# ----------------------------------------------------------------------------
# Passerelle MQTT d'un patient
# ----------------------------------------------------------------------------
class Gateway:
    """Client MQTT d'un patient : LWT, reconnexion auto, tampon hors-ligne."""

    def __init__(self, args, pid, events):
        import paho.mqtt.client as mqtt
        self.mqtt, self.args, self.pid, self.events = mqtt, args, pid, events
        self.gw_id = f"GW-{pid}"
        self.status_topic = f"{ROOT}/devices/{self.gw_id}/status"
        self.connected = False
        self.was_down = False
        self.buffer = deque(maxlen=MAX_BUFFER)
        self.client = None
        self.start()

    def start(self):
        mqtt, a = self.mqtt, self.args
        self.client = mqtt.Client(mqtt.CallbackAPIVersion.VERSION2,
                                  client_id=f"{self.gw_id}-{os.getpid()}")
        self.client.username_pw_set(a.user, a.password)
        if a.tls:
            self.client.tls_set(ca_certs=a.cafile)
        # Testament : si la passerelle disparait brutalement, le broker publie "offline"
        self.client.will_set(self.status_topic, "offline", qos=1, retain=True)
        self.client.reconnect_delay_set(min_delay=1, max_delay=5)
        self.client.on_connect = self._on_connect
        self.client.on_disconnect = self._on_disconnect
        self.client.connect_async(a.host, a.port, keepalive=10)
        self.client.loop_start()

    def _on_connect(self, client, userdata, flags, reason_code, properties):
        if reason_code.is_failure:
            print(f"[mqtt:{self.pid}] connexion refusee : {reason_code}")
            self.connected = False
            return
        self.connected = True
        client.publish(self.status_topic, "online", qos=1, retain=True)
        if self.was_down:
            self.events.log(self.pid, "broker_reconnected", f"{len(self.buffer)} en attente")
            self.was_down = False

    def _on_disconnect(self, client, userdata, disconnect_flags, reason_code, properties):
        self.connected = False
        if not self.was_down:
            self.was_down = True
            self.events.log(self.pid, "broker_disconnected", str(reason_code))

    def send(self, topic, payload):
        data = json.dumps(payload)
        if self.connected:
            self.flush()
            self.client.publish(topic, data, qos=1)
        else:
            self.buffer.append((topic, data))

    def flush(self):
        while self.buffer and self.connected:
            topic, data = self.buffer.popleft()
            self.client.publish(topic, data, qos=1)

    def crash(self):
        """Coupure brutale (sans DISCONNECT) : le broker doit publier le LWT."""
        self.client.loop_stop()
        sock = self.client.socket()
        if sock:
            sock.close()
        self.connected = False

    def close(self):
        """Arret propre : statut offline explicite (le LWT ne part pas sur un DISCONNECT normal)."""
        self.flush()
        if self.connected:
            info = self.client.publish(self.status_topic, "offline", qos=1, retain=True)
            info.wait_for_publish(timeout=2)
            self.client.disconnect()
        self.client.loop_stop()


# ----------------------------------------------------------------------------
# Programme principal
# ----------------------------------------------------------------------------
def parse_args():
    ap = argparse.ArgumentParser(description="Simulateur CTSMC (diabete) -> MQTT")
    ap.add_argument("--scenario", choices=SCENARIOS.keys(), default="normal")
    ap.add_argument("--patients", type=int, default=3, choices=[1, 2, 3])
    ap.add_argument("--interval", type=float, default=5.0, help="secondes entre deux pas")
    ap.add_argument("--duration", type=float, default=240.0, help="duree totale en s (0 = infini)")
    ap.add_argument("--anomaly-start", type=float, default=40.0)
    ap.add_argument("--anomaly-duration", type=float, default=90.0, help="duree du plateau (s)")
    ap.add_argument("--ramp", type=float, default=20.0, help="duree des rampes (s)")
    ap.add_argument("--loss-rate", type=float, default=0.3)
    ap.add_argument("--glitch-rate", type=float, default=0.4)
    ap.add_argument("--host", default="localhost")
    ap.add_argument("--port", type=int, default=1883)
    ap.add_argument("--user", default="simulator")
    ap.add_argument("--password", default=os.environ.get("MQTT_PASSWORD"),
                    help="sinon variable MQTT_PASSWORD, sinon saisie masquee")
    ap.add_argument("--tls", action="store_true")
    ap.add_argument("--cafile", default=None, help="certificat de l'autorite (avec --tls)")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--events-file",
                    default=str(Path(__file__).resolve().parent.parent / "results" / "events.csv"))
    ap.add_argument("--dry-run", action="store_true", help="affiche sans broker")
    return ap.parse_args()


def main():
    args = parse_args()
    if not args.dry_run and not args.password:
        args.password = getpass.getpass(f"Mot de passe MQTT pour '{args.user}' : ")

    rng = random.Random(args.seed)
    sc = SCENARIOS[args.scenario]
    run_id = f"{args.scenario}-{datetime.now().strftime('%Y%m%d-%H%M%S')}"
    events = EventLog(args.events_file, run_id)
    events.log("-", "run_start",
               f"scenario={args.scenario} seed={args.seed} interval={args.interval}")

    pids = list(PATIENTS)[:args.patients]
    models = {pid: PatientModel(pid, PATIENTS[pid], rng) for pid in pids}
    gws = {} if args.dry_run else {pid: Gateway(args, pid, events) for pid in pids}
    seq = {}
    target_pid = sc.get("patient") if sc else None
    if target_pid and target_pid not in models:
        print(f"[!] le scenario vise {target_pid} : utilise --patients 3")

    anomaly_end = args.anomaly_start + args.anomaly_duration
    started = ended = False
    t0 = time.monotonic()
    k = 0
    print(f"[sim] scenario={args.scenario} patients={pids} pas={args.interval}s (Ctrl+C pour arreter)")

    try:
        while True:
            t = time.monotonic() - t0
            if args.duration and t >= args.duration:
                break

            # Debut / fin d'anomalie (verite terrain)
            if sc and sc["kind"] != "manual":
                label = target_pid or "ALL"
                if not started and t >= args.anomaly_start:
                    events.log(label, "episode_start", f"type={args.scenario}")
                    print(f"[sim] >>> debut anomalie ({args.scenario}) a t={t:.0f}s")
                    started = True
                    if sc["kind"] == "crash" and target_pid in gws:
                        gws[target_pid].crash()
                        events.log(target_pid, "gateway_crash", "coupure brutale (LWT attendu)")
                if started and not ended and t >= anomaly_end:
                    events.log(label, "episode_end", f"type={args.scenario}")
                    print(f"[sim] <<< fin anomalie a t={t:.0f}s")
                    ended = True
                    if sc["kind"] == "crash" and target_pid in gws:
                        gws[target_pid].start()
                        events.log(target_pid, "gateway_restart", "")
            in_episode = started and not ended
            factor = ramp_factor(t, args.anomaly_start, args.anomaly_duration, args.ramp)

            line = []
            for pid in pids:
                targets = {}
                consent = True
                crashed = glitch = loss = False
                if sc and sc.get("patient") in (None, pid):
                    if sc["kind"] == "shift":
                        targets = sc["targets"]
                    elif sc["kind"] == "crash" and in_episode:
                        crashed = True
                    elif sc["kind"] == "glitch" and in_episode:
                        glitch = True
                    elif sc["kind"] == "consent" and in_episode:
                        consent = False
                    elif sc["kind"] == "loss" and in_episode:
                        loss = True

                values = models[pid].step(targets, factor if targets else 0.0, t)
                g_mgdl = values["glucose"]
                line.append(f"{pid}: {g_mgdl} mg/dL ({g_mgdl / 100:.2f} g/L)")

                for group, g in GROUPS.items():
                    if k % g["every"] != 0:
                        continue
                    key = (pid, group)
                    seq[key] = seq.get(key, 0) + 1   # le numero avance meme si le message est perdu
                    if crashed:
                        continue
                    if loss and rng.random() < args.loss_rate:
                        events.log(pid, "msg_dropped", f"{group} seq={seq[key]}")
                        continue
                    vals = dict(values)
                    if glitch and group == "glycemie" and rng.random() < args.glitch_rate:
                        vals["glucose"] = rng.choice([0, 1500, -10])
                        events.log(pid, "glitch_injected", f"value={vals['glucose']} seq={seq[key]}")
                    msg = build_message(pid, group, vals, consent)
                    msg["seq"] = seq[key]
                    topic = f"{ROOT}/patients/{pid}/{group}"
                    if args.dry_run:
                        print(topic, json.dumps(msg))
                    else:
                        gws[pid].send(topic, msg)

            print(f"t={t:5.0f}s  " + " | ".join(line))
            k += 1
            time.sleep(args.interval)
    except KeyboardInterrupt:
        print("\n[sim] arret demande")
    finally:
        for gw in gws.values():
            gw.close()
        events.log("-", "run_end", f"duration={time.monotonic() - t0:.0f}s")
        events.close()
        print("[sim] termine")


if __name__ == "__main__":
    sys.exit(main())