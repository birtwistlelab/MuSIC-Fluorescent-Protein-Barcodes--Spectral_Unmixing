import traceback
import importlib.util

step_files = [
    "my_module",
    "Step1_DataExtractionForReference_pos_pR_fp",
    "Step2_Profile_Reference",
    "Step3_DataExtractionForProbeUnmixing_pR_fp_singlets",
    "Step4_1_EntirePositivePopulationGating_FI_vs_frequency_hist",
    "Step4_2_SeparateTrainingAndTestingData",
    "Step5_UnmixTransfectedTrainingCells",
    "Step6_DetermineUnmixingThresholdsOnTrainingCells_FigS4",
    "Step7_UnmixingTransfectedTestingCells",
    "Step8_IdentifyEachCellInTestingDict_Fig3C",
    "Step9_1_RangeIntenistyCutoff",
    "Step9_2_CalculationAtEachCutoff",
    "Step9_3_Misclassified_mTFP1Cells_Fig3D",
    "Step9_4_OverallF1Score_vs_RemainingCells_Fig3E",
    "Step9_5_F1score_Comparision_Between_FPs_Fig3F",
    "Step10_DataExtractionForUnmixing_pos_pMuSIC",
    "Step11_IntensityCutoff_pMuSICs",
    "Step12_pMuSICs_Calculate_Correct_Fraction",
    "Step13_1_UnmixTransfectedTrainingCells_Remove_mTFP1",
    "Step13_2_DetermineTheThreshold_Remove_mTFP1",
    "Step13_3_Remove_mTFP1",
    "Step14_pMuSICs_Calculate_TPR_remove_mTFP1",
    "Step15_Reference_Spectrum_for_FigS2",
    "Step16_Troubleshooting_FigS6",
    "Step17_SpectraOverlapsBetweenCellLines_FigS3",
    "Step18_DifficultBarcodesBetweenCellLines_FigS7"
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