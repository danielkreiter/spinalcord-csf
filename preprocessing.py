# Import dependencies
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from peewee import *
from database.orm import *
from playhouse.shortcuts import model_to_dict, dict_to_model
from datetime import datetime, timedelta, date
from collections import OrderedDict
import seaborn as sns
mysql_db.close()

biological_effect = {"ALE" : 5 * 365, "OCR" : 0.5*365, "RTX" : 0.5*365, "CLA" : 96 * 7, "NAT" : 60, "OFA" : 0.5 * 365, "MIT" : 0.5 * 365} # Biological effectiveness of DMTs in days

hdmt = ["NAT", "OCR", "MIT", "RTX", "CLA", "ALE", "OFA"] # High efficacy-therapies
idmt = ["DMF", "DRF", "SIP", "FIN", "PON", "OZA"] # Intermediate efficacy-therapies
ldmt = ["GLA", "MTX", "AZA", "IVIG", "IFN", "TER"] # Low efficacy-therapies

# Cluster different inteferon variants
def replace_name(name):
    if "IFN" in name or "PEG" in name:
        return "IFN"
    else:
        return name

# Returns highest used DMT category within follow-up timeframe
def highest_category(patient, first_date, last_date):
    highest = 0
    highest_name = ""
    treatments = patient.treatments.where((reTreatment.name != "STUD") & (reTreatment.name != "PLAS") & (reTreatment.name != "PRED"))\
            .order_by(reTreatment.start_date.asc())
    for index, treatment in enumerate(treatments):

        if treatment.end_date == None:
            treatment.end_date = last_date

        if (treatment.name in ldmt or idmt) and abs(treatment.end_date - treatment.start_date).days < 90:
            continue

        if treatment.name in biological_effect.keys():
            if index+1 < len(treatments):
                treatment.end_date = min([treatment.end_date + timedelta(days=biological_effect[treatment.name]), treatments[index+1].start_date, last_date])
            else:
                treatment.end_date = min(treatment.end_date + timedelta(days=biological_effect[treatment.name]), last_date)
        if treatment.end_date > first_date and treatment.start_date < last_date:
            if treatment.name in hdmt:
                if highest < 3:
                    highest_name = treatment.name
                highest = np.max([3, highest])
            elif treatment.name in idmt:
                highest = np.max([2, highest])
                if highest < 2:
                    highest_name = treatment.name
            elif replace_name(treatment.name) in ldmt:
                highest = np.max([1, highest])
                if highest < 1:
                    highest_name = treatment.name
    return highest

def edss_around_date(patient, date):
    edss_col = []
    edss_col_date = []
    for edss in patient.edss:
        if abs(edss.date - date).days < 180:
            edss_col.append(edss.edss)
            edss_col_date.append(edss.date)
    if len(edss_col) >0:
        edss = min(edss_col)
        edss_date = edss_col_date[edss_col.index(edss)]
    else:
        edss = np.nan
        edss_date = np.nan
    return edss, edss_date

def relapses(patient, date_from, years):
    relapse_count = 0
    for relapse in patient.relapses:
        if (date_from - relapse.date).days < years*365 and (date_from - relapse.date).days > 0:
            relapse_count += 1
    return relapse_count

def brainmri_around_date(patient, date):
    baseline_brain_lesions = None
    baseline_brainmri_active = None
    brain_mri_query = patient.rebrainmris.where((reBrainMRI.date <= (date + timedelta(days=183))) & (reBrainMRI.date >= (date - timedelta(days=183))))
    interval = timedelta(365)
    match = None
    if brain_mri_query.count() > 0:
        baseline_brainmri_active = False
        for brain_i, brain_mri_match in enumerate(brain_mri_query.order_by(reBrainMRI.date)):
            n_interval = abs(date - brain_mri_match.date)
            if brain_mri_match.disease_activity:
                baseline_brainmri_active = True
            if n_interval < interval:
                match = brain_mri_match
                interval = n_interval
        baseline_brain_lesions = match.total_lesions
    return (baseline_brain_lesions, baseline_brainmri_active)

def pred_before(patient, date, interval):
    treatments = patient.treatments.where(reTreatment.name == "PRED").order_by(reTreatment.start_date.asc())
    for index, treatment in enumerate(treatments):
        if treatment.end_date:
            if (treatment.end_date - date).days > -interval and (treatment.end_date - date).days < (treatment.end_date - treatment.start_date).days:
                print("CSF {}, PRED END {}, {}".format(date, treatment.end_date, (date - treatment.end_date).days))
                return True
    return False

# Parameters of hyperbolic function (Reiber, 1998)
IgG_ab = 0.93
IgG_b2 = 6e-6
IgG_c = 1.7e-3

IgM_ab = 0.67
IgM_b2 = 120e-6
IgM_c = 7.1e-3

# In/exclusion parameters
max_timedelta_firstscmri = 366 # Days


# Create event tables
patients = rePatient.select().where(rePatient.post_exclusion == 0)

sc_date = []
sc_interval = []
sc_disease_activity = []
sc_newt2 = []
sc_newt1gd = []
sc_n_lesions = []
sc_last_mri_date = []

brain_date = []
brain_interval = []
brain_disease_activity = []
brain_newt2 = []
brain_newt1gd = []
brain_last_mri_date = []

pt_highest_category = []
pt_highest_category_brain = []
pt_edss = []
pt_relapses1y = []
pt_relapses2y = []
pt_brainmri_lesions = []
pt_brainmri_active = []
pt_age_at_baseline = []
pt_time_since_cis = []

eventtable_patient = []
eventtable_start = []
eventtable_end = []
eventtable_outcome = []
eventtable_thoracic = []

eventtable_brain_patient = []
eventtable_brain_start = []
eventtable_brain_end = []
eventtable_brain_outcome = []

count_i = 0

biological_effect = {"ALE" : 5*365, "OCR" : 0.5*365, "RTX" : 0.5*365, "CLA" : 96*7, 'OFA': 96*7, 'NAT' : 60}
levels = ['C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7', 'Th1', 'Th2', 'Th3', 'Th4', 'Th5', 'Th6', 'Th7', 'Th8', 'Th9', 'Th10', 'Th11', 'Th12']

# Still needs correction for first thoracic

for patient in patients:
    csf = patient.csf.get_or_none()
    if csf:
        if pred_before(patient, csf.date, 14):
            print("PRED BEFORE, EXCLUSION", patient.id)
            csf = False
    interval = None
    date = None
    disease_activity = None
    newt2 = None
    newt1gd = None
    n_lesions = None
    last_mri_date = None
    first_thoracic = None
    last_thoracic = None

    interval_brain = None
    date_brain = None
    disease_activity_brain = None
    newt2_brain = None
    newt1gd_brain = None
    last_brain_mri_date = None

    # If CSF data present, add data on first MRI to table
    if csf:
        spinal_mris_csf = patient.remris.where((reMRI.cervical == 1) & ((reMRI.date < csf.date + timedelta(days=max_timedelta_firstscmri)) & (reMRI.date > csf.date - timedelta(days=max_timedelta_firstscmri)))).order_by(reMRI.date)
        spinal_mris = patient.remris.where((reMRI.cervical == 1)).order_by(reMRI.date)
        brain_mris_csf = patient.rebrainmris.where(((reBrainMRI.date < csf.date + timedelta(days=max_timedelta_firstscmri)) & (reBrainMRI.date > csf.date - timedelta(days=max_timedelta_firstscmri)))).order_by(reBrainMRI.date)
        brain_mris = patient.rebrainmris.order_by(reBrainMRI.date)
        if spinal_mris_csf:

            # Spinal MRI
            df_spinal = pd.DataFrame(list(spinal_mris_csf.dicts()))
            df_spinal["interval"] = abs(df_spinal["date"] - csf.date)
            #if not spinal_mris[index].date > (patient.date_onset_cis - timedelta(days=30))
            first_mri = df_spinal[df_spinal["interval"] == df_spinal["interval"].min()]
            key = list(first_mri["interval"].keys())[0]
            interval = first_mri["interval"][key].days
            date = first_mri["date"][key]
            disease_activity = first_mri["disease_activity"][key]
            newt2 = first_mri["new_t2"][key]
            newt1gd = first_mri["new_gado"][key]
            n_lesions = reMRI.get_by_id(first_mri["id"][key]).lesions.count()
            th_mri_lesion = 0
            for lesion in reMRI.get_by_id(first_mri["id"][key]).lesions:
                if levels.index(lesion.level.split("-")[0]) >= levels.index("Th4"):
                    th_mri_lesion += 1
            n_lesions -= th_mri_lesion
            start = 0
            first_mri_index = 9999
            follow_up = 0

            if spinal_mris.count() > 1:
                # If more than 1 spinal cord MRI present
                for index, spinal_mri in enumerate(spinal_mris):
                    not_new_thoracic_lesions = 0
                    event_thoracic = 1
                    if spinal_mri.thoracic == 1 and not first_thoracic:
                        first_thoracic = index
                    elif spinal_mri.thoracic and first_thoracic:
                        if index > first_thoracic:
                            last_thoracic = index
                    else:
                        event_thoracic = 0

                    if spinal_mris[index].id == first_mri["id"][key]:
                        first_mri_index = index
                    elif index > first_mri_index:
                        eventtable_patient.append(patient.id)
                        eventtable_start.append(start)
                        end = abs(spinal_mris[index].date - first_mri["date"][key]).days
                        eventtable_end.append(end)
                        start = end
                        event_lesions = reMRI.get_by_id(spinal_mris[index].id).lesions.count()

                        # If this is a first spinal cord MRI with thoracic coverage, substract 'new' thoracic lesions
                        if spinal_mri.thoracic and index == first_thoracic:
                            for lesion in reMRI.get_by_id(spinal_mris[index].id).lesions:
                                if levels.index(lesion.level.split("-")[0]) >= levels.index("Th4") and lesion.new_t1_gado == 0:
                                    not_new_thoracic_lesions += 1
                        event_lesions -= not_new_thoracic_lesions
                        eventtable_thoracic.append(event_thoracic)
                        if event_lesions == 0:
                            eventtable_outcome.append(0)
                        else:
                            # eventtable_outcome.append(1)
                            if event_lesions > 1:
                                for i in range(event_lesions-1):
                                    eventtable_patient.append(patient.id)
                                    eventtable_start.append(start)
                                    end = start+0.1
                                    eventtable_end.append(end)
                                    eventtable_outcome.append(1)
                                    start = end

                        last_mri_date = spinal_mri.date
                        if index:
                            follow_up = 1
                            if index and not (spinal_mris[index].date or spinal_mris[index].date == ""):
            if follow_up == 1:
                count_i += 1

            # Brain MRI
            start_brain = 0
            first_mri_brain_index = 9999
            if len(list(brain_mris_csf.dicts())) > 0:
                df_brain = pd.DataFrame(list(brain_mris_csf.dicts()))
                if "date" not in df_brain.keys():
                    print(list(brain_mris_csf.dicts()))
                df_brain["interval"] = abs(df_brain["date"] - csf.date)

                first_mri_brain = df_brain[df_brain["interval"] == df_brain["interval"].min()]
                key_brain = list(first_mri_brain["interval"].keys())[0]
                interval_brain = first_mri_brain["interval"][key_brain].days
                date_brain = first_mri_brain["date"][key_brain]
                disease_activity_brain = first_mri_brain["disease_activity"][key_brain]
                newt2_brain = first_mri_brain["new_t2"][key_brain]
                newt1gd_brain = first_mri_brain["new_gado"][key_brain]

                if brain_mris.count() > 1:
                    for index, brain_mri in enumerate(brain_mris):
                        if brain_mris[index].id == first_mri_brain["id"][key_brain]:
                            first_mri_brain_index = index
                        elif index > first_mri_brain_index:
                            eventtable_brain_patient.append(patient.id)
                            eventtable_brain_start.append(start_brain)
                            end_brain = abs(brain_mris[index].date - first_mri_brain["date"][key_brain]).days
                            eventtable_brain_end.append(end_brain)
                            start_brain = end_brain
                            eventtable_brain_outcome.append(reBrainMRI.get_by_id(brain_mris[index].id).disease_activity)
                        last_brain_mri_date = brain_mri.date


    sc_date.append(date)
    sc_interval.append(interval)
    sc_disease_activity.append(disease_activity)
    sc_newt2.append(newt2)
    sc_newt1gd.append(newt1gd)
    sc_n_lesions.append(n_lesions)
    sc_last_mri_date.append(last_mri_date)

    #brain_date.append(date_brain)
    #brain_interval.append(interval_brain)
    #brain_disease_activity.append(disease_activity_brain)
    #brain_newt2.append(newt2_brain)
    #brain_newt1gd.append(newt1gd_brain)
    #brain_last_mri_date.append(last_brain_mri_date)

    # Add highest treatment category
    if last_mri_date:
        hc = highest_category(patient, date, last_mri_date)
    else:
        hc = None
    pt_highest_category.append(hc)

    if last_brain_mri_date:
        hc_brain = highest_category(patient, date, last_brain_mri_date)
    else:
        hc_brain = None
    pt_highest_category_brain.append(hc_brain)

    if csf:
        pt_edss.append(edss_around_date(patient, csf.date)[0])
        pt_relapses1y.append(relapses(patient, csf.date, 1))
        pt_relapses2y.append(relapses(patient, csf.date, 2))
        brain_mri_lesions, brain_mri_active = brainmri_around_date(patient, csf.date)
        pt_brainmri_lesions.append(brain_mri_lesions)
        pt_brainmri_active.append(brain_mri_active)
        age_at_baseline = int(patient.age - abs((datetime(2021, 5, 1).date() - csf.date)).days/365.25)
        pt_age_at_baseline.append(age_at_baseline)
        time_since_cis = (patient.date_onset_cis - csf.date).days/365.25
        pt_time_since_cis.append(time_since_cis)
    else:
        pt_edss.append(None)
        pt_relapses1y.append(None)
        pt_relapses2y.append(None)
        pt_brainmri_lesions.append(None)
        pt_brainmri_active.append(None)
        pt_age_at_baseline.append(None)
        pt_time_since_cis.append(None)

df_patients = pd.DataFrame(list(patients.dicts()))
df_patients["FIRST_SCMRI_DATE"] = sc_date
df_patients["FIRST_SCMRI_INTERVAL"] = sc_interval
df_patients["FIRST_SCMRI_ACTIVITY"] = sc_disease_activity
df_patients["FIRST_SCMRI_NEWT2LESIONS"] = sc_newt2
df_patients["FIRST_SCMRI_NEWGDLESIONS"] = sc_newt1gd
df_patients["FIRST_SCMRI_NEWLESIONCOUNT"] = sc_n_lesions
df_patients["LAST_SCMRI_DATE"] = sc_last_mri_date
df_patients["HIGHEST_CATEGORY"] = pt_highest_category
df_patients["EDSS"] = pt_edss
df_patients["RELAPSES_1Y"] = pt_relapses1y
df_patients["RELAPSES_2Y"] = pt_relapses2y
df_patients["BRAINMRI_LESIONS"] = pt_brainmri_lesions
df_patients["BRAINMRI_ACTIVE"] = pt_brainmri_active
df_patients["AGE_AT_BASELINE"] = pt_age_at_baseline
df_patients["TIME_SINCE_CIS"] = pt_time_since_cis

df_eventtable = pd.DataFrame()
df_eventtable["patient_id"] = eventtable_patient
df_eventtable["start"] = eventtable_start
df_eventtable["end"] = eventtable_end
df_eventtable["outcome"] = eventtable_outcome
df_eventtable["thoracic"] = eventtable_thoracic
df_eventtable.to_csv("csf_eventtable.csv")

#df_eventtable = pd.DataFrame()
#df_eventtable["patient_id"] = eventtable_brain_patient
#df_eventtable["start"] = eventtable_brain_start
#df_eventtable["end"] = eventtable_brain_end
#df_eventtable["outcome"] = eventtable_brain_outcome
#df_eventtable.to_csv("csf_eventtable_brain.csv")

# Patient data
query = reCSF.select()
df_csf = pd.DataFrame(list(query.dicts())) # Converteer resultaat van query naar dataframe
df_csf['date'] = pd.to_datetime(df_csf['date'], format="%Y-%m-%d")
n_total = len(df_patients.index)
n_null = len(df_patients[df_patients["FIRST_SCMRI_DATE"].isnull()].index)
n_both = len(df_patients[~df_patients["FIRST_SCMRI_DATE"].isnull() & ~df_patients["LAST_SCMRI_DATE"].isnull()].index)
print("{} without first SC MRI, {} with one, {} with also follow-up MRI".format(n_null, n_total-n_null, n_both))
df_patients_csf = df_patients.merge(df_csf, left_on="id", right_on = "patient_id")
df_patients_csf.head(10)

n_igg_complete = len(df_patients_csf[~np.isnan(df_patients_csf["igg"]) & ~np.isnan(df_patients_csf["alb"]) & ~np.isnan(df_patients_csf["index"]) & ~df_patients_csf["FIRST_SCMRI_DATE"].isnull() & (~np.isnan(df_patients_csf["alb_serum"]) | ~np.isnan(df_patients_csf["igg_serum"]))])
n_igm_complete = len(df_patients_csf[~np.isnan(df_patients_csf["igm"]) & ~np.isnan(df_patients_csf["alb"]) & ~np.isnan(df_patients_csf["index"]) & ~df_patients_csf["FIRST_SCMRI_DATE"].isnull() & (~np.isnan(df_patients_csf["alb_serum"]) | ~np.isnan(df_patients_csf["igg_serum"]))])
n_complete = len(df_patients_csf[~np.isnan(df_patients_csf["igm"]) & ~np.isnan(df_patients_csf["igg"]) & df_patients_csf["ocb"] & ~np.isnan(df_patients_csf["alb"]) & ~np.isnan(df_patients_csf["index"]) & ~df_patients_csf["FIRST_SCMRI_DATE"].isnull() & (~np.isnan(df_patients_csf["alb_serum"]) | ~np.isnan(df_patients_csf["igg_serum"]))])

print("{} complete CSF cases, {} IgG complete cases, {} IgM complete cases".format(n_complete, n_igg_complete, n_complete))
df_patients_csf_wbaseline = df_patients_csf[~df_patients_csf["FIRST_SCMRI_DATE"].isnull()]
df_patients_csf_wbaseline["Q_Alb"] = df_patients_csf_wbaseline["alb"] / df_patients_csf_wbaseline["alb_serum"]
df_patients_csf_wbaseline["Q_IgM"] = df_patients_csf_wbaseline["igm"] / df_patients_csf_wbaseline["igm_serum"]
df_patients_csf_wbaseline["Q_IgG"] = df_patients_csf_wbaseline["igg"] / df_patients_csf_wbaseline["igg_serum"]
df_patients_csf_wbaseline = df_patients_csf_wbaseline.reset_index()


def Qlim_IgG(x):
    return IgG_ab * np.sqrt(x**2+IgG_b2) - IgG_c
def Qlim_IgM(x):
    return IgM_ab * np.sqrt(x**2+IgM_b2) - IgM_c

def IF_IgG(QAlb, QIgG):
    Qlim = IgG_ab * np.sqrt(QAlb**2+IgG_b2) - IgG_c
    return (1-Qlim/QIgG)*100
def IF_IgM(QAlb, QIgM):
    Qlim = IgM_ab * np.sqrt(QAlb**2+IgM_b2) - IgM_c
    return (1-Qlim/QIgM)*100

x = np.linspace(1e-3, 25e-3, 100)
plt.plot(x, Qlim_IgM(x))
for i, val in enumerate(df_patients_csf_wbaseline["Q_Alb"]):
    plt.plot(df_patients_csf_wbaseline["Q_Alb"][i], df_patients_csf_wbaseline["Q_IgM"][i], marker="o", markersize=5)

#for i, txt in enumerate(df_patients_csf_wbaseline["IF_IgM"]):
#    plt.annotate(round(txt,1), (df_patients_csf_wbaseline["Q_Alb"][i], df_patients_csf_wbaseline["Q_IgM"][i]))
plt.show()

df_patients_csf_wbaseline["IF_IgG"] = IF_IgG(df_patients_csf_wbaseline["Q_Alb"], df_patients_csf_wbaseline["Q_IgG"])
df_patients_csf_wbaseline["IF_IgM"] = IF_IgM(df_patients_csf_wbaseline["Q_Alb"], df_patients_csf_wbaseline["Q_IgM"])
df_patients_csf_wbaseline['IgG_positive'] = np.where(df_patients_csf_wbaseline['IF_IgG'] > 0, 1, 0)
df_patients_csf_wbaseline.loc[(df_patients_csf_wbaseline['IgG_positive'] == 0) & (df_patients_csf_wbaseline['IF_IgG'].isnull()), 'IgG_positive'] = np.nan
df_patients_csf_wbaseline['IgM_positive'] = np.where(df_patients_csf_wbaseline['IF_IgM'] > 0, 1, 0)
df_patients_csf_wbaseline.loc[(df_patients_csf_wbaseline['IgM_positive'] == 0) & (df_patients_csf_wbaseline['IF_IgM'].isnull()), 'IgM_positive'] = np.nan


import re
def parse_ocb(ocb):
    m = re.findall(r"(\d+)-(\d+)", str(ocb))
    if len(m) > 0:
        ocb = ocb.split("-")[0]
    if not ocb:
        ocb = np.nan
    elif ocb in ['-', '0', '1']:
        ocb = 0
    elif ocb in ['?', '']:
        ocb = np.nan
    elif ocb in ['+']:
        ocb = 1
    elif ocb in ['>10', '10+', 10]:
        ocb = 2
    elif int(ocb) >= 10:
        ocb = 2
    elif int(ocb) < 10:
        ocb = 1
    else:
        ocb = np.nan
    return ocb

df_patients_csf_wbaseline["ocb_parsed"] = df_patients_csf_wbaseline["ocb"].apply(parse_ocb)
df_patients_csf_wbaseline.to_csv("csf_patients.csv")
