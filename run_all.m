function run_all()
%RUN_ALL  Run every numerical experiment for the paper and save all figures.
%
%   USAGE (MATLAB Online):  cd into this 'matlab' folder, then type:  run_all
%
%   DO NOT run the helper files directly at the >> prompt:
%       deg_sphere, assembled_field, reg_field, exportfig, local_ghat
%   They are FUNCTIONS that take arguments and are called only by s1..s6.
%   Typing their name alone gives a harmless 'Not enough input arguments'.
%
%   Scripts and what they check (paper cross-reference in parentheses):
%     S1  wedge-product magnitude identity            (eq. wedge)          PROVED
%     S2  planar CMN index = regularized index         (Thm cmn)           PROVED
%     S3  Gauss-map degree normalization on S^2         (eq. W)             PROVED
%     S4  candidate index (n=3) vs regularization       (eq. J/index)       EVIDENCE
%     S5  candidate index vs delta and delta->0 limit   (Conj. deltaindep)  EVIDENCE
%     S6  proximity condition as sufficient nonvanish   (eq. proximity)     PROVED
%
%   S1-S3, S6 verify PROVED statements. S4-S5 are NUMERICAL EVIDENCE only;
%   integrality of the candidate index for n>=3 remains an open problem.
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
