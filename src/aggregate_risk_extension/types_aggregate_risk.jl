
#######################################################
#              AGGREGATE RISK TYPES                   #
#######################################################

####################################################
# Technology
####################################################

## Mutable struct of Cobb_Douglas

mutable struct CobbDouglasTechnologyAR{R}
    const α::R 
    A::R 
    const δ::R 
    const μ::R 
end

#######################################################
# Household Workspace in periods 1 and 2
#######################################################

struct HouseholdWorkspace_T12{T1, S}
    lst_u_c
    lst_a
    lst_x
    lst_v
    lst_x_tmp
    lst_v_tmp
    lst_R_a_tmp

    u_c_tmp :: T1
    theta_star :: T1
    theta_pol :: T1
    a_pol_H :: T1 # This is a_3[cih_2,z_2] when A_2 is high. Note it is a function of the state [cih_2]
    a_pol_L :: T1 # This is a_3[cih_2,z_2] when A_2 is low. Note it is a function of the state [cih_2]
    k_pol :: T1  # This is k_2[a_1,z_1]. The standard policy function a_2[a_1,z_1] is stored in path[1,1].a_pol (it is the same in path[1,2].a_pol)
    b_pol :: T1  # This is b_2[a_1,z_1]. The standard policy function a_2[a_1,z_1] is stored in path[1,1].a_pol (it is the same in path[1,2].a_pol)
    resids_mat :: T1
    lower_index_H :: S
    lower_weight_H :: T1
    lower_index_L :: S
    lower_weight_L :: T1
    pdf_H :: T1 # This is pdf in the space [cih_2,z_2] if A is high
    pdf_L :: T1 # This is pdf in the space [cih_2,z_2] if A is low
end


function HouseholdWorkspace_T12(h::Household)  
    (; lst_u_c, lst_a, lst_x, lst_v, lst_x_tmp, lst_v_tmp, lst_R_a_tmp, u_c_tmp, theta_star, theta_pol, a_pol_H, a_pol_L, k_pol, b_pol, resids_mat,lower_index_H, lower_weight_H,
    lower_index_L,lower_weight_L, pdf_H, pdf_L) = _generate_base_workspace_matrices_T12(h)
    return HouseholdWorkspace_T12(lst_u_c, lst_a, lst_x, lst_v, lst_x_tmp, lst_v_tmp, lst_R_a_tmp, u_c_tmp, theta_star, theta_pol, a_pol_H, a_pol_L, k_pol, b_pol, resids_mat, lower_index_H, lower_weight_H,
    lower_index_L,lower_weight_L, pdf_H, pdf_L)
end 


function _generate_base_workspace_matrices_T12(h::Household)
    lst_u_c = Array{Any}(undef, (length(h.z_grid),2)) # n_z idiosyncratic states, and 2 aggregate productivity states: high and low.
    lst_a = similar(lst_u_c) # To be used in the forward iteration
    lst_x = similar(lst_u_c) 
    lst_v = similar(lst_u_c) 
    lst_x_tmp = Array{Any}(undef, (1,2))
    lst_v_tmp = similar(lst_x_tmp)
    lst_R_a_tmp = similar(lst_x_tmp)

    u_c_tmp = Array{eltype(h.a_grid)}(undef, length(h.a_grid), length(h.z_grid))
    theta_star = zeros(length(h.a_grid), length(h.z_grid))
    theta_pol = similar(theta_star)
    a_pol_H = similar(theta_star)
    a_pol_L = similar(theta_star)
    k_pol = similar(theta_star)
    b_pol = similar(theta_star)
    resids_mat = zeros(length(h.a_grid), length(h.z_grid),2)
    lower_index_H = similar(theta_star, Int)
    lower_weight_H = similar(theta_star)
    lower_index_L = similar(theta_star, Int)
    lower_weight_L = similar(theta_star)
    pdf_H = similar(theta_star)
    pdf_L = similar(theta_star)

    return (; lst_u_c, lst_a, lst_x, lst_v, lst_x_tmp, lst_v_tmp, lst_R_a_tmp, 
    u_c_tmp, theta_star, theta_pol, a_pol_H, a_pol_L, k_pol, b_pol, resids_mat, lower_index_H, lower_weight_H,
    lower_index_L,lower_weight_L, pdf_H, pdf_L)
end