import os
import traceback

step_files = [
    "my_module.py",
    "Step1_DataExtractionForReference(pos_pR_fp).py",
    "Step2_Profile_Reference.py",
    "Step3_DataExtractionForProbeUnmixing(pR_fp_singlets).py",
    "Step4_EntirePositivePopulationGating(FI_vs_frequency_hist).py",
    "Step5_UnmixPositiveCells(pR_fp_singlets).py",
    "Step6_Create18ROCCurves.py",
    "Step7_DetermineTheOptimal_ThresholdsForEachProbe.py",
    "Step8_1_Condition1_Original.py",
    "Step8_2_Identify_Each_Cell.py",
    "Step8_3_Fig.3D.py",
    "Step9.Condition2_Intensity_Cutoff.py",
    "Step10_1_Condition3_Remove_mTFP1.py",
    "Step10_2_determine_the_threshold_without_mTFP1.py",
    "Step10_3_FPR.py",
    "Step11_Condition4_IntensityCutoff_mTFP1Removal.py",
    "Step12_FPR_summary_of_4_conditions.py",
    "Step13_DataExtractionForUnmixing(pos_pMuSICs).py",
    "Step14_pMuSICs_condition4_IntensityCutoff_mTFP1Removal.py",
    "Step15_pMuSICs_condition4_Calculate_TPR.py",
    "Step16_Trouble_shooting_for_fig.S3.py",
    "Step17_Reference_Spectrum_for_fig.S2.py"
]


def run_all_steps():
    global_namespace = globals()
    local_namespace = {}
    global_namespace['np'] = __import__('numpy')
    for step in step_files:
        try:
            print(f"Running {step}...")
            with open(step) as file:
                exec(file.read(), global_namespace, local_namespace)
            print(f"{step} completed.\n")
        except Exception as e:
            print(f"Error while running {step}: {e}")
            traceback.print_exc()
            break


if __name__ == "__main__":
    run_all_steps()
