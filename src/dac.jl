function invcm_to_nm(shift_inv_cm::Vector{Float64}; laser_nm=532.0)
    laser_inv_cm = 1.0 ./ (laser_nm .* 1.0e-7)
    x_inv_cm = laser_inv_cm .- shift_inv_cm
    x_nm = 1.0 ./ (x_inv_cm .* 1.0e-7)
    return x_nm
end 
function nm_to_invcm(x::Vector{Float64}; laser_nm=532.0)
    x_inv_cm = 1.0./(x .*1.0e-7)
    laser_inv_cm = 1.0./(laser_nm*1.0e-7)
    shift_inv_cm = laser_inv_cm.-x_inv_cm
    return shift_inv_cm
end 