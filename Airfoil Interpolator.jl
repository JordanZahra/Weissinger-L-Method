
# Input Data
E395_Upper=[0 0 0.00004 0.00081 0.00029 0.00225 0.00073 0.00378 0.00175 0.00623 0.00323 0.00891 0.01088 0.01833 0.02281 0.0285 0.03891 0.03899 0.05906 0.04948 0.08309 0.0597 0.11082 0.06943 0.142 0.07845 0.17636 0.0866 0.21358 0.09372 0.25331 0.09967 0.29516 0.1043 0.33876 0.10748 0.38375 0.10908 0.42976 0.10904 0.47645 0.10732 0.52346 0.10395 0.57045 0.09897 0.6171 0.09248 0.66317 0.0847 0.70831 0.07601 0.752 0.06678 0.79373 0.05733 0.83295 0.04797 0.86911 0.03894 0.90168 0.03048 0.93013 0.02277 0.95395 0.01584 0.97298 0.0096 0.98734 0.0044 0.99669 0.00108 1 0]
E395_Lower=[0 0 0.00002 -0.00052 0.00011 -0.00112 0.0003 -0.00171 0.00058 -0.0023 0.00093 -0.0029 0.00183 -0.00411 0.00295 -0.00534 0.00504 -0.0072 0.01154 -0.01138 0.02505 -0.01685 0.04336 -0.02134 0.06641 -0.02468 0.09408 -0.02674 0.12626 -0.0275 0.16274 -0.02696 0.2033 -0.0252 0.24766 -0.02234 0.29548 -0.01859 0.34631 -0.01421 0.3996 -0.00949 0.45472 -0.00466 0.51099 0.00003 0.5677 0.00436 0.62409 0.00812 0.67936 0.01112 0.73272 0.01324 0.78334 0.01437 0.83041 0.01446 0.87314 0.01356 0.91076 0.01174 0.94258 0.00914 0.96783 0.00601 0.98584 0.00296 0.99649 0.00078 1 0]
E377_Upper=[0 0 0.00136 0.00482 0.00602 0.01244 0.01428 0.02102 0.02625 0.03013 0.04196 0.03945 0.06134 0.04871 0.08428 0.05767 0.11062 0.06611 0.14014 0.07381 0.17259 0.08059 0.20766 0.08625 0.24503 0.09054 0.28446 0.09316 0.32584 0.09397 0.36901 0.09292 0.41379 0.09006 0.45998 0.08553 0.50741 0.07959 0.55567 0.07272 0.60416 0.06527 0.65229 0.05751 0.69944 0.04966 0.74501 0.04193 0.78841 0.03449 0.82909 0.02751 0.8665 0.02112 0.90015 0.01545 0.92957 0.0106 0.95434 0.00664 0.97407 0.00363 0.9884 0.00157 0.99709 0.00039 1 0]
E377_Lower=[0 0 0.00324 -0.00406 0.01141 -0.00443 0.02442 -0.00204 0.0429 0.00363 0.06755 0.01263 0.09928 0.02461 0.13876 0.03866 0.18637 0.05342 0.24199 0.06713 0.30484 0.07754 0.37233 0.08211 0.44016 0.08077 0.50591 0.0758 0.56901 0.06879 0.62899 0.06056 0.68544 0.05179 0.738 0.04292 0.78634 0.03436 0.83016 0.02639 0.86923 0.01926 0.90333 0.01316 0.93231 0.00821 0.95607 0.00456 0.97474 0.00226 0.98849 0.001 0.99707 0.00028 1 0]
E171_Upper=[0 0 0.00297 0.00613 0.01138 0.01314 0.02461 0.02033 0.04247 0.02739 0.06481 0.03412 0.09142 0.04034 0.12207 0.04593 0.15645 0.05078 0.19422 0.05482 0.23501 0.05796 0.27839 0.06013 0.32389 0.06127 0.37102 0.06123 0.41943 0.05984 0.46887 0.05707 0.51907 0.05302 0.56975 0.04795 0.62048 0.04228 0.67062 0.03638 0.7195 0.03052 0.76641 0.02493 0.81068 0.0198 0.8516 0.01525 0.88853 0.01134 0.92086 0.00808 0.94806 0.00536 0.96985 0.00302 0.98613 0.0012 0.99644 0.00023 1 0]
E171_Lower=[0 0 0.00297 -0.00613 0.01138 -0.01314 0.02461 -0.02033 0.04247 -0.02739 0.06481 -0.03412 0.09142 -0.04034 0.12207 -0.04593 0.15645 -0.05078 0.19422 -0.05482 0.23501 -0.05796 0.27839 -0.06013 0.32389 -0.06127 0.37102 -0.06123 0.41943 -0.05984 0.46887 -0.05707 0.51907 -0.05302 0.56975 -0.04795 0.62048 -0.04228 0.67062 -0.03638 0.7195 -0.03052 0.76641 -0.02493 0.81068 -0.0198 0.8516 -0.01525 0.88853 -0.01134 0.92086 -0.00808 0.94806 -0.00536 0.96985 -0.00302 0.98613 -0.0012 0.99644 -0.00023 1 0]

# Variables
StepSize_x=(2/3)*0.01 # SET TO 0.01 for 100 cuts, (2/3)*0.01 for 151 cuts
StepSize_y=(1/3)*0.01 # SET TO 0.005 for 201 cuts, 0.004 for 251 cuts,
Added_Thickness=0.0025 # Added thickness to the E377 Airfoil (this is half of the percentage chord added thickness)
Dihedral_Angle_Deg=6 # Quarter chord dihedral angle in degrees (note that dihedral is positive and anhedral is negative)
Smoothing_Addition_Index=0.11 # SHOULD BE 0.11 OR GREATER TO ALIGN WITH FUSELAGE TERM APPROACHING ZERO BUT BETTER RESUTLS AT 0.05 TO 0.06 ISH
Dihedral_Step=0.001
J3=100000
Voting_Integer=42 # (28 for StepSize_y=0.005, 37 for StepSize_y=0.0025 (both for StepSize_x=0.01) (42 for StepSize_y=0.004 for StepSize_x=(2/3)*0.02)

# Data manipulation block - E395
Length_E395_Upper=Int64(length(E395_Upper)/2)
Length_E395_Lower=Int64(length(E395_Lower)/2)
E395_Upper_xz=reshape(E395_Upper,2,Length_E395_Upper)
E395_Upper_x=E395_Upper_xz[1,:]
E395_Upper_z=E395_Upper_xz[2,:]
E395_Lower_xz=reshape(E395_Lower,2,Length_E395_Lower)
E395_Lower_x=E395_Lower_xz[1,:]
E395_Lower_z=E395_Lower_xz[2,:]

# Data manipulation block - E377
Length_E377_Upper=Int64(length(E377_Upper)/2)
Length_E377_Lower=Int64(length(E377_Lower)/2)
E377_Upper_xz=reshape(E377_Upper,2,Length_E377_Upper)
E377_Upper_x=E377_Upper_xz[1,:]
E377_Upper_z=E377_Upper_xz[2,:]
E377_Lower_xz=reshape(E377_Lower,2,Length_E377_Lower)
E377_Lower_x=E377_Lower_xz[1,:]
E377_Lower_z=E377_Lower_xz[2,:]

# Data manipulation block - E171
Length_E171_Upper=Int64(length(E171_Upper)/2)
Length_E171_Lower=Int64(length(E171_Lower)/2)
E171_Upper_xz=reshape(E171_Upper,2,Length_E171_Upper)
E171_Upper_x=E171_Upper_xz[1,:]
E171_Upper_z=E171_Upper_xz[2,:]
E171_Lower_xz=reshape(E171_Lower,2,Length_E171_Lower)
E171_Lower_x=E171_Lower_xz[1,:]
E171_Lower_z=E171_Lower_xz[2,:]

# Linearization block - E395
Gradient_E395_Upper=zeros(Length_E395_Upper-1)
Offset_E395_Upper=zeros(Length_E395_Upper-1)
for i=1:Length_E395_Upper-1
Gradient_E395_Upper[i]=(E395_Upper_z[1+i]-E395_Upper_z[i])/(E395_Upper_x[1+i]-E395_Upper_x[i])
Offset_E395_Upper[i]=E395_Upper_z[i]-Gradient_E395_Upper[i]*E395_Upper_x[i]
end
Gradient_E395_Lower=zeros(Length_E395_Lower-1)
Offset_E395_Lower=zeros(Length_E395_Lower-1)
for i=1:Length_E395_Lower-1
Gradient_E395_Lower[i]=(E395_Lower_z[1+i]-E395_Lower_z[i])/(E395_Lower_x[1+i]-E395_Lower_x[i])
Offset_E395_Lower[i]=E395_Lower_z[i]-Gradient_E395_Lower[i]*E395_Lower_x[i]
end

# Linearization block - E377
Gradient_E377_Upper=zeros(Length_E377_Upper-1)
Offset_E377_Upper=zeros(Length_E377_Upper-1)
for i=1:Length_E377_Upper-1
Gradient_E377_Upper[i]=(E377_Upper_z[1+i]-E377_Upper_z[i])/(E377_Upper_x[1+i]-E377_Upper_x[i])
Offset_E377_Upper[i]=E377_Upper_z[i]-Gradient_E377_Upper[i]*E377_Upper_x[i]
end
Gradient_E377_Lower=zeros(Length_E377_Lower-1)
Offset_E377_Lower=zeros(Length_E377_Lower-1)
for i=1:Length_E377_Lower-1
Gradient_E377_Lower[i]=(E377_Lower_z[1+i]-E377_Lower_z[i])/(E377_Lower_x[1+i]-E377_Lower_x[i])
Offset_E377_Lower[i]=E377_Lower_z[i]-Gradient_E377_Lower[i]*E377_Lower_x[i]
end

# Linearization block - E171
Gradient_E171_Upper=zeros(Length_E171_Upper-1)
Offset_E171_Upper=zeros(Length_E171_Upper-1)
for i=1:Length_E171_Upper-1
Gradient_E171_Upper[i]=(E171_Upper_z[1+i]-E171_Upper_z[i])/(E171_Upper_x[1+i]-E171_Upper_x[i])
Offset_E171_Upper[i]=E171_Upper_z[i]-Gradient_E171_Upper[i]*E171_Upper_x[i]
end
Gradient_E171_Lower=zeros(Length_E171_Lower-1)
Offset_E171_Lower=zeros(Length_E171_Lower-1)
for i=1:Length_E171_Lower-1
Gradient_E171_Lower[i]=(E171_Lower_z[1+i]-E171_Lower_z[i])/(E171_Lower_x[1+i]-E171_Lower_x[i])
Offset_E171_Lower[i]=E171_Lower_z[i]-Gradient_E171_Lower[i]*E171_Lower_x[i]
end

# Airfoil template construction and replication - E395
E395_Upper_x_Template=0:StepSize_x:1
Length_E395_Upper_Template=Int64(length(E395_Upper_x_Template))
E395_Upper_z_Template=zeros(Length_E395_Upper_Template)
for i=1:Length_E395_Upper_Template
for j=1:Length_E395_Upper-1
if (E395_Upper_x_Template[i]>=E395_Upper_x[j])&&(E395_Upper_x_Template[i]<E395_Upper_x[j+1])
E395_Upper_z_Template[i]=Gradient_E395_Upper[j]*E395_Upper_x_Template[i]+Offset_E395_Upper[j]
end
end
end
E395_Lower_x_Template=0:StepSize_x:1
Length_E395_Lower_Template=Int64(length(E395_Lower_x_Template))
E395_Lower_z_Template=zeros(Length_E395_Lower_Template)
for i=1:Length_E395_Lower_Template
for j=1:Length_E395_Lower-1
if (E395_Lower_x_Template[i]>=E395_Lower_x[j])&&(E395_Lower_x_Template[i]<E395_Lower_x[j+1])
E395_Lower_z_Template[i]=Gradient_E395_Lower[j]*E395_Lower_x_Template[i]+Offset_E395_Lower[j]
end
end
end

# Airfoil template construction and replication - E377
E377_Upper_x_Template=0:StepSize_x:1
Length_E377_Upper_Template=Int64(length(E377_Upper_x_Template))
E377_Upper_z_Template=zeros(Length_E377_Upper_Template)
for i=1:Length_E377_Upper_Template
for j=1:Length_E377_Upper-1
if (E377_Upper_x_Template[i]>=E377_Upper_x[j])&&(E377_Upper_x_Template[i]<E377_Upper_x[j+1])
E377_Upper_z_Template[i]=Gradient_E377_Upper[j]*E377_Upper_x_Template[i]+Offset_E377_Upper[j]
end
end
end
E377_Lower_x_Template=0:StepSize_x:1
Length_E377_Lower_Template=Int64(length(E377_Lower_x_Template))
E377_Lower_z_Template=zeros(Length_E377_Lower_Template)
for i=1:Length_E377_Lower_Template
for j=1:Length_E377_Lower-1
if (E377_Lower_x_Template[i]>=E377_Lower_x[j])&&(E377_Lower_x_Template[i]<E377_Lower_x[j+1])
E377_Lower_z_Template[i]=Gradient_E377_Lower[j]*E377_Lower_x_Template[i]+Offset_E377_Lower[j]
end
end
end

# Airfoil template construction and replication - E171
E171_Upper_x_Template=0:StepSize_x:1
Length_E171_Upper_Template=Int64(length(E171_Upper_x_Template))
E171_Upper_z_Template=zeros(Length_E171_Upper_Template)
for i=1:Length_E171_Upper_Template
for j=1:Length_E171_Upper-1
if (E171_Upper_x_Template[i]>=E171_Upper_x[j])&&(E171_Upper_x_Template[i]<E171_Upper_x[j+1])
E171_Upper_z_Template[i]=Gradient_E171_Upper[j]*E171_Upper_x_Template[i]+Offset_E171_Upper[j]
end
end
end
E171_Lower_x_Template=0:StepSize_x:1
Length_E171_Lower_Template=Int64(length(E171_Lower_x_Template))
E171_Lower_z_Template=zeros(Length_E171_Lower_Template)
for i=1:Length_E171_Lower_Template
for j=1:Length_E171_Lower-1
if (E171_Lower_x_Template[i]>=E171_Lower_x[j])&&(E171_Lower_x_Template[i]<E171_Lower_x[j+1])
E171_Lower_z_Template[i]=Gradient_E171_Lower[j]*E171_Lower_x_Template[i]+Offset_E171_Lower[j]
end
end
end

# Repair of the E377 Airfoil
Difference_UpperLowerE377=E377_Upper_z_Template-E377_Lower_z_Template
Difference_UpperLowerE377[findin(E377_Upper_x_Template,E377_Upper_x_Template[E377_Upper_x_Template.>=0.9])]=0.01
Difference_UpperLowerE377[findin(E377_Upper_x_Template,E377_Upper_x_Template[E377_Upper_x_Template.<=0.1])]=0.01
Minimum_Difference_UpperLowerE377=minimum(Difference_UpperLowerE377)
Location_Minimum_Difference_UpperLowerE377=findin(Difference_UpperLowerE377,Minimum_Difference_UpperLowerE377)
Location_Minimum_Difference_UpperLowerE377=Location_Minimum_Difference_UpperLowerE377[1]
for i=Location_Minimum_Difference_UpperLowerE377:Length_E171_Lower_Template
if (E377_Upper_z_Template[i]-E377_Lower_z_Template[i])>=Minimum_Difference_UpperLowerE377
E377_Lower_z_Template[i]=E377_Upper_z_Template[i]
end
end
E377_Upper_z_Template=E377_Upper_z_Template+Added_Thickness
E377_Lower_z_Template=E377_Lower_z_Template-Added_Thickness
E377_Upper_z_Template[1]=0
E377_Lower_z_Template[1]=0
E377_Upper_z_Template[Length_E377_Upper_Template]=0
E377_Lower_z_Template[Length_E377_Lower_Template]=0

# Interpolation on the output planform
m_Symm=Int64((m+1)/2)
Gradient_c_v=zeros(m_Symm-1)
Offset_c_v=zeros(m_Symm-1)
Gradient_c_v_FUSELAGE=zeros(m_Symm-1)
Offset_c_v_FUSELAGE=zeros(m_Symm-1)
Gradient_x_AC=zeros(m_Symm-1)
Offset_x_AC=zeros(m_Symm-1)
Gradient_Geometric_Twist=zeros(m_Symm-1)
Offset_Geometric_Twist=zeros(m_Symm-1)
for i=1:m_Symm-1
Gradient_c_v[i]=(c_v_WASHOUT[1+i]-c_v_WASHOUT[i])/(Nu_v_WASHOUT[1+i]-Nu_v_WASHOUT[i])
Offset_c_v[i]=c_v_WASHOUT[i]-Gradient_c_v[i]*Nu_v_WASHOUT[i]
Gradient_c_v_FUSELAGE[i]=(c_v_FUSELAGE[1+i]-c_v_FUSELAGE[i])/(Nu_v_WASHOUT[1+i]-Nu_v_WASHOUT[i])
Offset_c_v_FUSELAGE[i]=c_v_FUSELAGE[i]-Gradient_c_v_FUSELAGE[i]*Nu_v_WASHOUT[i]
Gradient_x_AC[i]=(x_AC[1+i]-x_AC[i])/(Nu_v_WASHOUT[1+i]-Nu_v_WASHOUT[i])
Offset_x_AC[i]=x_AC[i]-Gradient_x_AC[i]*Nu_v_WASHOUT[i]
Gradient_Geometric_Twist[i]=(Geometric_Twist_Distribution_Deg_WASHOUT[1+i]-Geometric_Twist_Distribution_Deg_WASHOUT[i])/(Nu_v_WASHOUT[1+i]-Nu_v_WASHOUT[i])
Offset_Geometric_Twist[i]=Geometric_Twist_Distribution_Deg_WASHOUT[i]-Gradient_Geometric_Twist[i]*Nu_v_WASHOUT[i]
end
Wingspan_y_Template=0:StepSize_y:1
Length_Wingspan_y_Template=Int64(length(Wingspan_y_Template))
c_v_TEMPLATE=zeros(Length_Wingspan_y_Template)
c_v_FUSELAGE_TEMPLATE=zeros(Length_Wingspan_y_Template)
x_AC_TEMPLATE=zeros(Length_Wingspan_y_Template)
Geomtric_Twist_TEMPLATE=zeros(Length_Wingspan_y_Template)
for i=1:Length_Wingspan_y_Template
for j=1:m_Symm-1
if (Wingspan_y_Template[i]>=Nu_v_WASHOUT[j+1])&&(Wingspan_y_Template[i]<Nu_v_WASHOUT[j])
c_v_TEMPLATE[i]=Gradient_c_v[j]*Wingspan_y_Template[i]+Offset_c_v[j]
c_v_FUSELAGE_TEMPLATE[i]=Gradient_c_v_FUSELAGE[j]*Wingspan_y_Template[i]+Offset_c_v_FUSELAGE[j]
x_AC_TEMPLATE[i]=Gradient_x_AC[j]*Wingspan_y_Template[i]+Offset_x_AC[j]
Geomtric_Twist_TEMPLATE[i]=Gradient_Geometric_Twist[j]*Wingspan_y_Template[i]+Offset_Geometric_Twist[j]
end
end
end
c_v_TEMPLATE[1]=c_v_WASHOUT[m_Symm]
c_v_FUSELAGE_TEMPLATE[1]=c_v_FUSELAGE[m_Symm]
x_AC_TEMPLATE[1]=x_AC[m_Symm]
Geomtric_Twist_TEMPLATE[1]=Geometric_Twist_Distribution_Deg_WASHOUT[m_Symm]
Geomtric_Twist_TEMPLATE[length(Geomtric_Twist_TEMPLATE)]=Geometric_Twist_Distribution_Deg_WASHOUT[1]

# So what I need is the quarter chord x-location (sweep), y-location (winspan coordinate), and z-location (dihedral)
# At each of these locations, I then calculate the blended airfoil based on a blending function (function of y-location only)
# Then I plot the airfoil at the x,y,z location by taking the applicable blended airfoil, and plotting it in 2D (on the x-z plane) 25% forward in the x-direction and 75% backwards in the x-direction

# Quarter chord y-location
FINAL_y_locations=Wingspan_y_Template*b/2

# Quarter chord x-location
Quarter_Chord_Sweep_Line=FINAL_y_locations.*tan(-Quarter_Chord_Sweep_Angle_Deg*pi/180)
Quarter_Chord_Sweep_Line=Quarter_Chord_Sweep_Line+x_AC_TEMPLATE-c_v_FUSELAGE_TEMPLATE*Fuselage_x_Offset_Index
Quarter_Chord_Sweep_Line=Quarter_Chord_Sweep_Line-Quarter_Chord_Sweep_Line[1]
LE_Line=Quarter_Chord_Sweep_Line+0.25*c_v_TEMPLATE
TE_Line=Quarter_Chord_Sweep_Line-0.75*c_v_TEMPLATE
FINAL_x_locations=LE_Line

# Quarter chord z-location
# The rest is calculated post belding, scaling and twisting as it is dependent on this to occur first
Dihedral_Line=Wingspan_y_Template.*tan(Dihedral_Angle_Deg*pi/180)
Dihedral_Line_ORIGINAL=Dihedral_Line
h_Function=zeros(Length_Wingspan_y_Template)
h_Function_Location=zeros(Length_Wingspan_y_Template)
h_Function_ADDITIONAL=zeros(Length_Wingspan_y_Template)
dh_dy_Function=zeros(Length_Wingspan_y_Template,Length_E377_Upper_Template)
dh_dy_Plotter=zeros(Length_Wingspan_y_Template-1)
Airfoil_Slice_Location_Plotter=zeros(Length_Wingspan_y_Template)
FINAL_z_locations=Dihedral_Line

# Blending function
# Blending function
# Blending function

# Finding point at which AlphaL0=-5 degrees
AlphaL0_Distribution_Deg=(AlphaL0_Tip_Deg-AlphaL0_Root_Deg)*Wingspan_y_Template+AlphaL0_Root_Deg
AlphaL0_Negative5_Location=(-5-AlphaL0_Root_Deg)./(AlphaL0_Tip_Deg-AlphaL0_Root_Deg)
Wingspan_y_Template_Finder=Wingspan_y_Template.'.'
Wingspan_y_Template_Finder[findin(Wingspan_y_Template_Finder,Wingspan_y_Template_Finder[Wingspan_y_Template_Finder.<AlphaL0_Negative5_Location])]=0
AlphaL0_Negative5_Index=findfirst(Wingspan_y_Template_Finder)
if (Wingspan_y_Template[AlphaL0_Negative5_Index]-AlphaL0_Negative5_Location)<=(Wingspan_y_Template[AlphaL0_Negative5_Index-1]-AlphaL0_Negative5_Location)
AlphaL0_Negative5_Index=AlphaL0_Negative5_Index
else
AlphaL0_Negative5_Index=AlphaL0_Negative5_Index-1
end

Airfoil_XZ_Matrix_Upper=zeros(Length_Wingspan_y_Template,Length_E377_Upper_Template)
Airfoil_XZ_Matrix_Lower=zeros(Length_Wingspan_y_Template,Length_E377_Lower_Template)
R1=(0-1)/(Wingspan_y_Template[AlphaL0_Negative5_Index]-0)*Wingspan_y_Template+1
R2=(1-0)/(1-Wingspan_y_Template[AlphaL0_Negative5_Index])*Wingspan_y_Template+(1-(1-0)/(1-Wingspan_y_Template[AlphaL0_Negative5_Index]))

# Blending E395 to E377 from AlphaL0_Negaitve5_Location to tip [LINEAR BLEND]
for i=1:Length_Wingspan_y_Template
for j=1:Length_E377_Upper_Template
if Wingspan_y_Template[i]<=Wingspan_y_Template[AlphaL0_Negative5_Index]
Airfoil_XZ_Matrix_Upper[i,j]=R1[i]*E395_Upper_z_Template[j]+(1-R1[i])*E377_Upper_z_Template[j]
Airfoil_XZ_Matrix_Lower[i,j]=R1[i]*E395_Lower_z_Template[j]+(1-R1[i])*E377_Lower_z_Template[j]
end
end
end

# Blending E377 to E171 from root to AlphaL0_Negaitve5_Location [NON-LINEAR BLEND]
for i=1:Length_Wingspan_y_Template
for j=1:Length_E377_Lower_Template
if Wingspan_y_Template[i]>=Wingspan_y_Template[AlphaL0_Negative5_Index]
Airfoil_XZ_Matrix_Upper[i,j]=(1-R2[i])*E377_Upper_z_Template[j]+R2[i]*E171_Upper_z_Template[j]
Airfoil_XZ_Matrix_Lower[i,j]=(1-R2[i])*E377_Lower_z_Template[j]+R2[i]*E171_Lower_z_Template[j]
end
end
end

# Applying blending, scaling and washout
# Applying blending, scaling and washout
# Applying blending, scaling and washout
Airfoil_Scaling_Factor=c_v_TEMPLATE
Geomtric_Twist_TEMPLATE=-Geomtric_Twist_TEMPLATE*pi/180
Rotation_Output_Upper_x=zeros(Length_Wingspan_y_Template,length(E395_Upper_x_Template))
Rotation_Output_Upper_z=zeros(Length_Wingspan_y_Template,length(E395_Upper_x_Template))
Rotation_Output_Lower_x=zeros(Length_Wingspan_y_Template,length(E395_Lower_x_Template))
Rotation_Output_Lower_z=zeros(Length_Wingspan_y_Template,length(E395_Lower_x_Template))
for i=1:Length_Wingspan_y_Template
Rotation_Center_x=0.25
Rotation_Center_z=0
Rotation_Center=repmat([Rotation_Center_x;Rotation_Center_z],1,length(E395_Upper_x_Template))
Rotation_Matrix=[cos(Geomtric_Twist_TEMPLATE[i]) -sin(Geomtric_Twist_TEMPLATE[i]); sin(Geomtric_Twist_TEMPLATE[i]) cos(Geomtric_Twist_TEMPLATE[i])]
Upper_XZ_Twisted=zeros(2,length(E395_Upper_x_Template))
Upper_XZ_Twisted[1,:]=Upper_XZ_Twisted[1,:]+E395_Upper_x_Template
Upper_XZ_Twisted[2,:]=Upper_XZ_Twisted[2,:]+Airfoil_XZ_Matrix_Upper[i,:]
Lower_XZ_Twisted=zeros(2,length(E395_Lower_x_Template))
Lower_XZ_Twisted[1,:]=Lower_XZ_Twisted[1,:]+E395_Lower_x_Template
Lower_XZ_Twisted[2,:]=Lower_XZ_Twisted[2,:]+Airfoil_XZ_Matrix_Lower[i,:]
Output_Upper_XZ_Twisted=Rotation_Matrix*(Upper_XZ_Twisted-Rotation_Center)+Rotation_Center
Output_Lower_XZ_Twisted=Rotation_Matrix*(Lower_XZ_Twisted-Rotation_Center)+Rotation_Center
Rotation_Output_Upper_x[i,:]=Output_Upper_XZ_Twisted[1,:]
Rotation_Output_Upper_z[i,:]=Output_Upper_XZ_Twisted[2,:]
Rotation_Output_Lower_x[i,:]=Output_Lower_XZ_Twisted[1,:]
Rotation_Output_Lower_z[i,:]=Output_Lower_XZ_Twisted[2,:]
end

# Final airfoils (scaled to millimetres not metres) for output and plotting
AIRFOIL_OUTPUT_x_Upper=zeros(Length_Wingspan_y_Template,Length_E377_Upper_Template)
AIRFOIL_OUTPUT_x_Lower=zeros(Length_Wingspan_y_Template,Length_E377_Lower_Template)
AIRFOIL_OUTPUT_y=zeros(Length_Wingspan_y_Template,Length_E377_Upper_Template)
AIRFOIL_OUTPUT_z_Upper=zeros(Length_Wingspan_y_Template,Length_E377_Upper_Template)
AIRFOIL_OUTPUT_z_Lower=zeros(Length_Wingspan_y_Template,Length_E377_Lower_Template)
for i=1:Length_Wingspan_y_Template
for j=1:length(E377_Upper_x_Template)
AIRFOIL_OUTPUT_x_Upper[i,j]=Rotation_Output_Upper_x[i,j].*Airfoil_Scaling_Factor[i]-FINAL_x_locations[i]
AIRFOIL_OUTPUT_x_Lower[i,j]=Rotation_Output_Lower_x[i,j].*Airfoil_Scaling_Factor[i]-FINAL_x_locations[i]
AIRFOIL_OUTPUT_y[i,j]=FINAL_y_locations[i]
AIRFOIL_OUTPUT_z_Upper[i,j]=Rotation_Output_Upper_z[i,j]*Airfoil_Scaling_Factor[i]+FINAL_z_locations[i]
AIRFOIL_OUTPUT_z_Lower[i,j]=Rotation_Output_Lower_z[i,j]*Airfoil_Scaling_Factor[i]+FINAL_z_locations[i]
end
end
AIRFOIL_OUTPUT_x_Upper=AIRFOIL_OUTPUT_x_Upper-AIRFOIL_OUTPUT_x_Upper[1,1]+0.1
AIRFOIL_OUTPUT_x_Lower=AIRFOIL_OUTPUT_x_Lower-AIRFOIL_OUTPUT_x_Lower[1,1]+0.1
AIRFOIL_OUTPUT_x_Upper=AIRFOIL_OUTPUT_x_Upper*1000
AIRFOIL_OUTPUT_x_Lower=AIRFOIL_OUTPUT_x_Lower*1000
AIRFOIL_OUTPUT_y=AIRFOIL_OUTPUT_y*1000
AIRFOIL_OUTPUT_z_Upper=AIRFOIL_OUTPUT_z_Upper*1000
AIRFOIL_OUTPUT_z_Lower=AIRFOIL_OUTPUT_z_Lower*1000
Quarter_Chord_Sweep_Line=Quarter_Chord_Sweep_Line*1000
TE_Line=TE_Line*1000
LE_Line=LE_Line*1000
FINAL_y_locations=FINAL_y_locations*1000
Dihedral_Line=Dihedral_Line*1000



# Continuation of Quarter chord z-location
for i=1:Length_Wingspan_y_Template
h_Function[i]=maximum(AIRFOIL_OUTPUT_z_Upper[i,:])
Airfoil_Slice=AIRFOIL_OUTPUT_z_Upper[i,:]
Airfoil_Slice_Location=findin(Airfoil_Slice,Airfoil_Slice[Airfoil_Slice.==h_Function[i]])
Airfoil_Slice_Location=Airfoil_Slice_Location[1]
h_Function_Location[i]=AIRFOIL_OUTPUT_x_Upper[i,Airfoil_Slice_Location]
Airfoil_Slice_Location_Plotter[i]=Airfoil_Slice_Location
end
for i=1:Length_Wingspan_y_Template-1
dh_dy_Function[i]=(h_Function[i+1]-h_Function[i])/(Wingspan_y_Template[i+1]-Wingspan_y_Template[i])
end
# START OF ITERATION LOOP
Adder_Integer=zeros(length(E377_Upper_x_Template))
Dihedral_Smoothing_y_Template=zeros(Length_Wingspan_y_Template)+Wingspan_y_Template
Dihedral_Smoothing_y_Template[findin(Dihedral_Smoothing_y_Template,Dihedral_Smoothing_y_Template[Dihedral_Smoothing_y_Template.>Smoothing_Addition_Index])]=0
Smoothing_Addition_IndexLocator=findlast(Dihedral_Smoothing_y_Template)
for k=1:J3
println("j3 = ",k," out of ",J3)
println(" ")


for i=1:Length_Wingspan_y_Template-1
for j=1:length(E377_Upper_x_Template)
dh_dy_Function[i,j]=(AIRFOIL_OUTPUT_z_Upper[i+1,j]-AIRFOIL_OUTPUT_z_Upper[i,j])/(Wingspan_y_Template[i+1]-Wingspan_y_Template[i])
dh_dy_Plotter[i]=(Wingspan_y_Template[i+1]+Wingspan_y_Template[i])/2
end
end


for i=1:Smoothing_Addition_IndexLocator
for j=1:length(E377_Upper_x_Template)
if (dh_dy_Function[i,j]<dh_dy_Function[Smoothing_Addition_IndexLocator+1,j])&&(dh_dy_Function[i,j]<dh_dy_Function[i+1,j])&&(j<=Airfoil_Slice_Location_Plotter[i])&&(AIRFOIL_OUTPUT_x_Upper[i,j]>AIRFOIL_OUTPUT_x_Upper[Smoothing_Addition_IndexLocator,1])
Adder_Integer[j]=1
end
end
if sum(Adder_Integer)>Voting_Integer # roughly 0.6*Airfoil_Slice_Location_Plotter[1]
h_Function_ADDITIONAL[i]=h_Function_ADDITIONAL[i]+Dihedral_Step
println("SUM ADDER INTEGER = ",sum(Adder_Integer),", i = ",i)
println(" ")
Adder_Integer=Adder_Integer.*0
end
end


for i=1:Length_Wingspan_y_Template
AIRFOIL_OUTPUT_z_Upper[i,:]=AIRFOIL_OUTPUT_z_Upper[i,:]-h_Function_ADDITIONAL[i]
AIRFOIL_OUTPUT_z_Lower[i,:]=AIRFOIL_OUTPUT_z_Lower[i,:]-h_Function_ADDITIONAL[i]
end
Dihedral_Line=Dihedral_Line-h_Function_ADDITIONAL
FINAL_z_locations=FINAL_z_locations-h_Function_ADDITIONAL
h_Function_ADDITIONAL=h_Function_ADDITIONAL.*0


end




using Plots
plotly()
plot()
xlims!(0,1)
ylims!(-1000,1000)
plot!(dh_dy_Plotter,dh_dy_Function)
gui()

using Plots
plotly()
plot()
plot!(Wingspan_y_Template,Dihedral_Line_ORIGINAL*1000)
plot!(Wingspan_y_Template,Dihedral_Line)
gui()

using Plots
plotly()
plot()
xlims!(-1600,1600)
ylims!(-1600,1600)
plot!(FINAL_y_locations,Quarter_Chord_Sweep_Line)
plot!(FINAL_y_locations,-h_Function_Location+325.11675)
plot!(FINAL_y_locations,LE_Line)
plot!(FINAL_y_locations,TE_Line)
plot!(-FINAL_y_locations,Quarter_Chord_Sweep_Line)
plot!(-FINAL_y_locations,-h_Function_Location+325.11675)
plot!(-FINAL_y_locations,LE_Line)
plot!(-FINAL_y_locations,TE_Line)
gui()

using Plots
plotly()
plot()
xlims!(-1600,1600)
ylims!(-1600,1600)
zlims!(-1600,1600)
for i=1:Length_Wingspan_y_Template
plot!(AIRFOIL_OUTPUT_x_Upper[i,:],AIRFOIL_OUTPUT_y[i,:],AIRFOIL_OUTPUT_z_Upper[i,:])
plot!(AIRFOIL_OUTPUT_x_Lower[i,:],AIRFOIL_OUTPUT_y[i,:],AIRFOIL_OUTPUT_z_Lower[i,:])
end
for i=1:Length_Wingspan_y_Template
plot!(AIRFOIL_OUTPUT_x_Upper[i,:],-AIRFOIL_OUTPUT_y[i,:],AIRFOIL_OUTPUT_z_Upper[i,:])
plot!(AIRFOIL_OUTPUT_x_Lower[i,:],-AIRFOIL_OUTPUT_y[i,:],AIRFOIL_OUTPUT_z_Lower[i,:])
end
plot!(-Quarter_Chord_Sweep_Line+325.11675,FINAL_y_locations,Dihedral_Line)
plot!(h_Function_Location,FINAL_y_locations,h_Function)
plot!(-LE_Line+325.11675,FINAL_y_locations,Dihedral_Line)
plot!(-TE_Line+325.11675,FINAL_y_locations,Dihedral_Line)
plot!(-Quarter_Chord_Sweep_Line+325.11675,-FINAL_y_locations,Dihedral_Line)
plot!(h_Function_Location,-FINAL_y_locations,h_Function)
plot!(-LE_Line+325.11675,-FINAL_y_locations,Dihedral_Line)
plot!(-TE_Line+325.11675,-FINAL_y_locations,Dihedral_Line)
gui()

using Plots
plotly()
plot()
plot!(E395_Upper_x_Template,E395_Upper_z_Template)
plot!(E395_Lower_x_Template,E395_Lower_z_Template)
plot!(E377_Upper_x_Template,E377_Upper_z_Template)
plot!(E377_Lower_x_Template,E377_Lower_z_Template)
plot!(E171_Upper_x_Template,E171_Upper_z_Template)
plot!(E171_Lower_x_Template,E171_Lower_z_Template)
gui()




#
