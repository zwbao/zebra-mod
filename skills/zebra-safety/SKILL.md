---
name: zebra-safety
description: Urgent warnings that must not wait for a research answer — when a rare-disease presentation needs emergency care today (metabolic decompensation, status epilepticus, adrenal crisis, cardiac or respiratory decline, stroke-like episodes), and the drugs and procedures that are dangerous in specific rare diseases (sodium-channel blockers in Dravet, anaesthesia in muscular dystrophy and mitochondrial disease, HLA-B*15:02 and carbamazepine, aminoglycosides in m.1555A>G). Read when a symptom is happening now, or for a treatment, drug, anaesthesia or procedure question — not for every symptom in a history. Triggers automatically from the other zebra skills; also on 急不急, 要不要去医院, 急诊, 能不能用这个药, 麻醉, 手术, contraindication, emergency, is this urgent.
---

# zebra-safety — what cannot wait, and what can harm

This skill says when to stop researching and seek care, and which specific hazards matter. The wording rules in sections 1–3 are for patients and families; a clinician is the care team and gets hazards, mechanism and acute management directly (section 4). The hazards below are established textbook contraindications: state them at once, cite a retrieved source when you have one (`disease_card`, `gene_card`, GeneReviews, `literature_search`), and when you have not checked this session, state them anyway labelled "to verify" — a hazard is never withheld because a lookup failed or has not run yet.

## 1. Stop and seek care now

If the person describes any of these as happening now or today, say so first, before any analysis, in plain words, and then continue only if they want to:

- **A seizure lasting more than 5 minutes, or seizures that repeat without recovery** between them (status epilepticus) — emergency care. First aid while help comes: on the side, nothing in the mouth except a rescue medicine the child's doctors prescribed for this (buccal or nasal midazolam goes into the cheek or nose — follow their written plan, never give a dose yourself), note the time.
- **Vomiting, lethargy, refusing food, rapid breathing or confusion in a child with (or suspected to have) a metabolic disorder**, especially during an infection, fasting or after a high-protein meal — possible metabolic decompensation; many inborn errors of metabolism have an emergency protocol, and hours matter. Ask whether they have an emergency letter or protocol from their metabolic team, and tell them to bring it.
- **Any sudden weakness, drooping face, trouble speaking, or a stroke-like episode** — emergency care (and in mitochondrial disease, MELAS-type episodes present this way).
- **Breathlessness at rest, a new inability to lie flat, morning headaches with daytime sleepiness** in a neuromuscular disease — possible respiratory failure; urgent assessment.
- **Fainting, a racing or very slow heartbeat, chest pain** in a disease with cardiac involvement (muscular dystrophies, mitochondrial disease, Fabry, connective-tissue disorders) — urgent.
- **Severe vomiting, low blood pressure, extreme weakness in adrenal insufficiency or a steroid-treated patient** — adrenal crisis; emergency care, and stress-dose steroid rules come from their own team.
- **A child on long-term steroids with an infection, injury or surgery** — they may need a stress dose; contact the treating team.
- **Any new severe pain, high fever with an implanted device or a central line, or a sudden change the family calls "not like the usual episodes"** — ask them to contact the team rather than wait for research.

Wording to a family: "This sounds like something to get seen for today — please contact your doctor or go to the emergency department; I can keep looking at the records meanwhile." To a family: never estimate how dangerous it is, never suggest a dose, never tell them to wait. To a clinician: say it is urgent, then give the acute management you would propose (first-line steps, emergency-protocol content, dose ranges as proposals) — they are the ones acting on it.

## 2. Disease-specific hazards worth raising

To a family, present each as a question for the care team, never as an instruction; to a clinician, as the hazard, its mechanism and what to do instead. Cite when checked; label "to verify" when not.

| Setting | Hazard to ask about |
|---|---|
| Dravet syndrome and other SCN1A loss-of-function epilepsies | Sodium-channel-blocking antiseizure medicines (carbamazepine, oxcarbazepine, phenytoin, lamotrigine) can make seizures worse; a history of worsening on one of these is itself a clue to the diagnosis |
| Any Han Chinese (and broader Asian) patient before carbamazepine or oxcarbazepine | HLA-B\*15:02 testing for the risk of severe skin reactions; pharmacogenomic testing is a question for the prescriber |
| Duchenne/Becker and other muscular dystrophies, myotonia, mitochondrial disease | Anaesthesia risk: depolarising muscle relaxants (succinylcholine) and volatile agents can cause rhabdomyolysis, hyperkalaemic cardiac arrest or malignant-hyperthermia-like reactions; the anaesthetist must know the diagnosis before any procedure, including dental sedation |
| Mitochondrial disease, m.1555A>G and related variants | Aminoglycoside antibiotics can cause profound, permanent hearing loss; also fasting and valproate are often avoided |
| Urea-cycle disorders and organic acidaemias | Valproate, prolonged fasting, high-protein loads and steroids can precipitate crises; emergency regimens exist (to a family: they come from the metabolic team; to a clinician: give the protocol's content as a proposal) |
| Long-QT and other arrhythmia syndromes | QT-prolonging drugs (some antibiotics, antiemetics, antipsychotics) — any new drug should be checked against a QT list |
| Marfan, Loeys-Dietz, vascular Ehlers-Danlos | Contact sports, isometric strain, fluoroquinolones; aortic imaging intervals |
| Osteogenesis imperfecta | Handling and positioning during procedures; bisphosphonate decisions belong to the specialist |
| Fabry, Pompe, other enzyme-replacement therapies | Infusion reactions and antibody formation — managed by the treating centre |
| Any gene-therapy or antisense trial | Eligibility often closes after certain treatments or antibody exposure; ask before starting anything that could exclude the child later |

## 3. How to say it

- To a family: one short paragraph, plain words, what to do and who to contact. No statistics, no mechanism, no reassurance you cannot support.
- To a clinician: the hazard, the mechanism, and the source.
- Record every hazard raised as a question in the case (`case_update` → `questions`) so it reaches the next visit.

## 4. What this skill never does

To a patient or family: dosing, stopping or starting a medicine, estimating prognosis for this person, or deciding whether an episode is "safe to watch at home". Those are the treating team's, and saying so is not a failure to help — it is the honest answer.

To a clinician the hazards above are not a ceiling: they get the full picture — mechanism, source, dose ranges and prognosis as proposals for their review (`zebra-care`). Urgent-first still holds for everyone.
