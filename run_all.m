function run_all()
    here = fileparts(mfilename('fullpath'));
    addpath(here);
    fprintf('\n########## Poincare-Hopf / Filippov  numerical suite ##########\n');
    s1_wedge_identity();
    s2_planar_index_cmn();
    s3_gauss_degree_S2();
    s4_candidate_index_n3();
    s5_delta_independence();
    s6_proximity_condition();
    fprintf('\n########## done ##########\n');
    fprintf('Copy the printed values above and send them back for comparison.\n');
end
