import os

step_files = [
    "src/my_module.py",
    "src/Step1_DataExtractionForReference(pos_pR_fp).py",
    "src/Step2_Profile_Reference.py",
    "src/Step3_DataExtractionForProbeUnmixing(pR_fp_singlets).py",
    "src/Step4_EntirePositivePopulationGating(FI_vs_frequency_hist).py",
    "src/Step5_UnmixPositiveCells(pR_fp_singlets).py",
    "src/Step6_Create18ROCCurves.py",
    "src/Step7_DetermineTheOptimal_ThresholdsForEachProbe.py",
    "src/Step8_1_Condition1_Original.py",
    "src/Step8_2_Identify_Each_Cell.py",
    "src/Step8_3_Fig.3D.py",
    "src/Step9.Condition2_Intensity_Cutoff.py",
    "src/Step10_1_Condition3_Remove_mTFP1.py",
    "src/Step10_2_determine_the_threshold_without_mTFP1.py",
    "src/Step10_3_FPR.py",
    "src/Step11_Condition4_IntensityCutoff_mTFP1Removal.py",
    "src/Step12_FPR_summary_of_4_conditions.py",
    "src/Step13_DataExtractionForUnmixing(pos_pMuSICs).py",
    "src/Step14_pMuSICs_condition4_IntensityCutoff_mTFP1Removal.py",
    "src/Step15_pMuSICs_condition4_Calculate_TPR.py",
    "src/Step16_Trouble_shooting_for_fig.S3.py",
    "src/Step17_Reference_Spectrum_for_fig.S2.py"
]


def run_all_steps():
    for step in step_files:
        try:
            print(f"Running {step}...")
            with open(step) as file:
                exec(file.read())
            print(f"{step} completed.\n")
        except Exception as e:
            print(f"Error while running {step}: {e}")
            break


if __name__ == "__main__":
    run_all_steps()
