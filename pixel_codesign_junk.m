%% transistor-level OTA design: action (one stage)

% M1: input pair (NMOS)

% vectors to keep track of results for comparison
W = zeros(num_lengths_N, 1);
WL = zeros(num_lengths_N, 1);
gmro = zeros(num_lengths_N, 1);
gmid = zeros(num_lengths_N, 1);

% try each length
for i = 1:num_lengths_N
    
    % column for gmid, gmro data
    col = 1 + i;
    
    % determine corresponding length
    L_charac = L_N(i);
    W_charac = L_charac * 2;
    
    for j = 1:points_per_length
        
        % go from SI to WI to get minimum required gmid
        point = points_per_length - j + 1;
        
        gmid(i) = gmid_n(point, col);
        
        if (gmid(i) > gmid_WI)
            % found required gmid
            
            % determine corresponding current density for this length
            ID_charac = gmid_n(point, 1);
            ID_density = ID_charac / W_charac;
            
            % determine actual ID required
            ID_req = gm_min / gmid(i);
            ID = 1e-6 * ceil(ID_req * 1e6);
            
            % determine actual W required
            W(i) = ID / ID_density;
            WL(i) = W(i) .* L_charac;
            
            % keep track of corresponding gmro
            gmro(i) = gmro_n(point, col);
            break
        end      
    end
end

% plot results for M1
figure
plot(L_N*1e6, W*1e6);
xlabel('L (\mum)');
ylabel('W (\mum)');
title('M1: W vs. L');
figure
plot(L_N*1e6, WL*1e12);
xlabel('L (\mum)');
ylabel('W ((\mum)^2)');
title('M1: WL vs. L');
figure
plot(L_N*1e6, gmro);
xlabel('L (\mum)');
ylabel('g_mr_o');
title('M1: g_mr_o vs. L');

% pick L that gives us required WL for parasitics/flicker
for i = 1:num_lengths_N
    WL_this = WL(i);
    
    if (WL_this > WL_min)
        % determine final dimensions
        W_M1 = W(i);
        L_M1 = L_N(i);
        % determine final small signal params
        gmid_M1 = gmid(i);
        gmro_M1 = gmro(i);
        break
    end
end

% compute some currents
ID_main = ID;
ID_load = ID_main / current_scale;
ID_both = ID_main + ID_load;

% compute M1 small signal params
gm_M1 = gmid_M1 * ID_main;
ro_M1 = gmro_M1 / gm_M1;

% beginning to compute the gain

% Gm assuming that not all current through load (pessimistic)
Gm_target = gm_M1 / gain_factor;
Rout_target = A_OL_target / Gm_target;

% assume Rout = gm9ro9ro11 (fact that this is not true is taken
% into account in gain_factor)
ro_11_target = Rout_target / gmro_9_target;

% M11: load device (NMOS)
% the only strong-inversion transistor (bc noise)
[W, gmid, gmro] = size_M(L_N, gmid_n, gmro_n, gmid_SI, 0, ID_load,...
    num_lengths_N, points_per_length);
% compute gmro target
gm_11_assumed = gmid_SI * ID_load;
gmro_11_target = gm_11_assumed * ro_11_target;
% pick L
[L_M11, W_M11, gmid_M11, gmro_M11] = pick_L(L_N, W, gmid, gmro,...
    gmro_11_target, num_lengths_N, W_min);
% compute M11 small signal params
gm_M11 = gmid_M11 * ID_load;
ro_M11 = gmro_M11 / gm_M11;

% recompute gmro target for M9 based on actual ro_M11
gmro_9_target = Rout_target / ro_M11;

% M9: load cascode device (NMOS)
[W, gmid, gmro] = size_M(L_N, gmid_n, gmro_n, gmid_WI, 1, ID_load,...
    num_lengths_N, points_per_length);
[L_M9, W_M9, gmid_M9, gmro_M9] = pick_L(L_N, W, gmid, gmro,...
    gmro_9_target, num_lengths_N, W_min);
% compute M9 small signal params
gm_M9 = gmid_M9 * ID_load;
ro_M9 = gmro_M9 / gm_M9;

% size the degen resistors for reasonable voltage swing
R_degen = V_R_degen / ID_both;

% determine (gmro) target for 3,5,7 based on relation that should be gg
r_down = gmro_M9 * ro_M11;
gmro_sq_target = r_down * gg_factor / R_degen;

gmro_357_target = sqrt(gmro_sq_target) * gmro_reduc;

% size 3,5,7 based on this target

% M3: telecascode device (NMOS)
[W, gmid, gmro] = size_M(L_N, gmid_n, gmro_n, gmid_WI, 1, ID_main,...
    num_lengths_N, points_per_length);
[L_M3, W_M3, gmid_M3, gmro_M3] = pick_L(L_N, W, gmid, gmro,...
    gmro_357_target, num_lengths_N, W_min);
% compute M3 small signal params
gm_M3 = gmid_M3 * ID_main;
ro_M3 = gmro_M3 / gm_M3;

% M5: current source device (PMOS)
[W, gmid, gmro] = size_M(L_P, gmid_p, gmro_p, gmid_WI, 1, ID_both,...
    num_lengths_P, points_per_length);
[L_M5, W_M5, gmid_M5, gmro_M5] = pick_L(L_P, W, gmid, gmro,...
    gmro_357_target, num_lengths_P, W_min);
% compute M5 small signal params
gm_M5 = gmid_M5 * ID_both;
ro_M5 = gmro_M5 / gm_M5;

% M7: foldcascode device (PMOS)
[W, gmid, gmro] = size_M(L_P, gmid_p, gmro_p, gmid_WI, 1, ID_load,...
    num_lengths_P, points_per_length);
[L_M7, W_M7, gmid_M7, gmro_M7] = pick_L(L_P, W, gmid, gmro,...
    gmro_357_target, num_lengths_P, W_min);
% compute M7 small signal params
gm_M7 = gmid_M7 * ID_load;
ro_M7 = gmro_M7 / gm_M7;

% calculate the actual gain

% impedance looking into M3 drain
r_inbranch = gmro_M3 * ro_M1;
% impedance looking into folding node
r_fold = 1 / (1/(gmro_M5 * R_degen) + 1/(r_inbranch));
% impedance looking up from output
r_up = r_fold * gmro_M7;

% output resistance
Rout_actual = 1/(1/r_down + 1/r_up);

% transconductance
Gs7 = (gm_M7 + 1/ro_M7) / (1 + gmro_M9 * ro_M11 / ro_M7);
Gd5 = 1 / (ro_M5 * (1 + R_degen / ro_M5 + gm_M5 * R_degen));
Gs3 = (gm_M3 + 1/ro_M3) / (1 + 1 / (ro_M3 * (Gs7 + Gd5)));

input_split = Gs3 * ro_M1 / (1 + Gs3 * ro_M1);
fold_split = Gs7 / (Gs7 + Gd5);

Gm_actual = gm_M1 * input_split * fold_split;

A_OL_actual = Rout_actual * Gm_actual;