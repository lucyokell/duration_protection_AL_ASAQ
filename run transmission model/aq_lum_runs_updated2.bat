#### RUN MODEL WHEN ASAQ PROPHYLAXIS IS LONG, AL SHORT.
set root=indiv_p_gam

set exe=%root%\check_indiv_p_gam_2.exe

set parms=default_parms.txt
set site_file=site_files\site_generic.txt
set output_root=output

REM ******general sim parms *********
set opts=output_header 1 output_type 0 num_runs 3 baseline_int 0 init_output 10 gamb_ss_c 1
set opts_seas=output_header 1 output_type 0 num_runs 3 baseline_int 0 init_output 10 gamb_ss_c 0.001

set opts_drugs= p_protect 0 drug0_eff 0.6 drug1_eff 0.95 drug1_rel_c .05094 drug2_eff 0.95 drug2_rel_c .102 drug1_dur_P 8.729359 add drug1_shape_P 93.49584 drug2_dur_P 17.80105 add drug2_shape_P 16.8051

set change_drug=drug_cov_0_0 0.8 drug_cov_1_0 0 change_drug 1 drug_cov_0_1 0 drug_cov_1_1 0.8
set recalculate=recalculate 2 prev_age0 2 prev_age1 10


REM ############## RUN MODEL ###################
REM ******** drug effects no seasonality **********
%exe% %root% %parms% %site_file% %output_root% run_lum_short_ns_5 %opts% %opts_drugs% %change_drug% %recalculate% prev 0.053 change_drug_time 0 final_run 6  num_people 600000
%exe% %root% %parms% %site_file% %output_root% run_aq_long_ns_5 %opts% %opts_drugs% %change_drug% %recalculate% prev 0.053 drug_cov_1_1 0  drug_cov_2_1 0.8  change_drug_time 0 final_run 6  num_people 600000

%exe% %root% %parms% %site_file% %output_root% run_lum_short_ns_15 %opts% %opts_drugs% %change_drug% %recalculate% prev 0.16   change_drug_time 0 final_run 6  num_people 600000
%exe% %root% %parms% %site_file% %output_root% run_aq_long_ns_15 %opts% %opts_drugs% %change_drug% %recalculate% prev 0.16   drug_cov_1_1 0  drug_cov_2_1 0.8  change_drug_time 0 final_run 6  num_people 600000

%exe% %root% %parms% %site_file% %output_root% run_lum_short_ns_50 %opts% %opts_drugs% %change_drug% %recalculate% prev 0.546   change_drug_time 0 final_run 6  num_people 600000
%exe% %root% %parms% %site_file% %output_root% run_aq_long_ns_50 %opts% %opts_drugs% %change_drug% %recalculate% prev 0.546   drug_cov_1_1 0  drug_cov_2_1 0.8  change_drug_time 0 final_run 6  num_people 600000

******** drug effect, with seasonality. Set prevalence to achieve the same EIR **********
%exe% %root% %parms% %site_file% %output_root% run_lum_short_s_5 %opts_seas% %opts_drugs% %change_drug% %recalculate% prev 0.0465 change_drug_time 0 final_run 6  num_people 600000
%exe% %root% %parms% %site_file% %output_root% run_aq_long_s_5 %opts_seas% %opts_drugs% %change_drug% %recalculate% prev 0.0465 drug_cov_1_1 0  drug_cov_2_1 0.8  change_drug_time 0 final_run 6   num_people 600000

%exe% %root% %parms% %site_file% %output_root% run_lum_short_s_15 %opts_seas% %opts_drugs% %change_drug% %recalculate% prev 0.1305 change_drug_time 0 final_run 6   num_people 600000
%exe% %root% %parms% %site_file% %output_root% run_aq_long_s_15 %opts_seas% %opts_drugs% %change_drug% %recalculate% prev 0.1305 drug_cov_1_1 0  drug_cov_2_1 0.8  change_drug_time 0 final_run 6   num_people 600000

%exe% %root% %parms% %site_file% %output_root% run_lum_short_s_50 %opts_seas% %opts_drugs% %change_drug% %recalculate% prev 0.379 change_drug_time 0 final_run 6  num_people 600000
%exe% %root% %parms% %site_file% %output_root% run_aq_long_s_50 %opts_seas% %opts_drugs% %change_drug% %recalculate% prev 0.379 drug_cov_1_1 0  drug_cov_2_1 0.8  change_drug_time 0 final_run 6  num_people 600000

pause