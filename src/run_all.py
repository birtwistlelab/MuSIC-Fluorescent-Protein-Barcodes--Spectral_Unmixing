import traceback
import importlib.util
import os

step_files = [
    "my_module",
    "Step1_DataExtractionForReference_pos_pR_fp",
    "Step2_Profile_Reference",
    "Step3_DataExtractionForProbeUnmixing_pR_fp_singlets",
    "Step4_EntirePositivePopulationGating_FI_vs_frequency_hist",
    "Step5_UnmixPositiveCells_pR_fp_singlets",
    "Step6_Create18ROCCurves",
    "Step7_DetermineTheOptimal_ThresholdsForEachProbe",
    "Step8_1_Condition1_Original",
    "Step8_2_Identify_Each_Cell",
    "Step8_3_Fig_3D",
    "Step9_Condition2_Intensity_Cutoff",
    "Step10_1_Condition3_Remove_mTFP1",
    "Step10_2_determine_the_threshold_without_mTFP1",
    "Step10_3_FPR",
    "Step11_Condition4_IntensityCutoff_mTFP1Removal",
    "Step12_FPR_summary_of_4_conditions",
    "Step13_DataExtractionForUnmixing_pos_pMuSICs",
    "Step14_pMuSICs_condition4_IntensityCutoff_mTFP1Removal",
    "Step15_pMuSICs_condition4_Calculate_TPR",
    "Step16_Trouble_shooting_for_fig_S3",
    "Step17_Reference_Spectrum_for_fig_S2"
]


def run_all_steps():
    for step in step_files:
        try:
            print(f"Running {step}...")
            module = importlib.import_module(step)
            print(f"{step} completed.\n")
            
        except Exception as e:
            print(f"Error while running {step}: {e}")
            traceback.print_exc()
            break


if __name__ == "__main__":
    run_all_steps()
