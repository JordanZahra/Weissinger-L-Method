
#
#
#
#
#
#
#
#
#
#

# CONSTRAINTS ARE AS FOLLOWS
# 1. Only calculates symmetrical spanloads about the root chord (this is important as there are NaNs and Infs in the code that will present themselves in the LStar_v_Mu computation if an asymmetric loading is attempted)
# 2. No discontinuities in washout are permissable
# 3. The quarter chord line must be straight across the semispan
# 4. The wing must be planar
# 5. The flow is assumed to be inviscid (inapppropriate for flows which are dominated by viscous forces or when the the wing begins to stall)
# 6. There is only one chordwise control point (inappropriate for low aspect ratio wings (less than approximately AR=6) or for highly swept wings (sweep angles greater than 45 degrees))
# 7. There is an inbuilt error associated with approximation using Fourier transformations that scale with 1/(m^2)

# PART 1 CONSTANTS
# Constants to be modified prior to running
# Note that there is a strong correlation with the AR and the sweep angle. Convergence will be difficult for high AR wings at extreme forward or aft sweep angles however, it is still possible to have extreme forward or aft sweep angles with low AR wings or vice versa.
# Note that there is difficulty in setting the AlphaL0_Distribution and Washout_Distribution before the grid independence study is complete as the m value will change, hence the detail required in the distributions will change. For this reason, the distributions must be updated below within the planform optimization code if they are not a constant
FreeStream_Velocity=10 # In metres per second (note that the faster the freestream velocity, the smaller the induced angle of attack)
Sigma_SpanRatio=4/3
m=599 # WARNING: MUST BE AN ODD NUMBER (so that one point is at the centreline of the distribution)
b=3 # Wingspan in Metres
Root_Chord=0.5 # Root chord in metres
Quarter_Chord_Sweep_Angle_Deg=30 # Quarter chord sweep angle in degrees (note that aft sweep is positive)
Mach_Number=FreeStream_Velocity/339.242732420313 # Note that a value of zero still provides valid results but does not account for the effects of compressibility
CL_RealWorld=2*pi # Set to 2*pi for theoretical value or change to the known 3D experiemental value
Design_AoA_Deg=5 # With reference to the root chord
Geometric_Twist_Distribution_Deg=0 # With reference to the root chord. Note that due to the way this method treats the airfoil as a flat plate, this is a modification to the AlphaL0_Distribution and as such a negative washout actually increases lift rather than decreasing it (this is accounted for in the Weissinger code by reversing the sign of the geometric twist when calculating Alpha)
AlphaL0_Root_Deg=-6.36363636363636
AlphaL0_Tip_Deg=0
Plotter_Step_Size=0.01

# PART 2 CONSTANTS
m_Limit_Computer=599 # SET THIS TO 599
m_Increment=2 # Must be divisible by 2
Accuracy_Increment=0.01 # 1% of the solution of maximum memory of the computer (assumed to be correct enough to use as a grid independent solution)

# PART 3 CONSTANTS
Planform_Adjust_Factor=0.0000001 # Accurate to 0.01mm
J1=15000 # Note that this needs to be larger for span loadings that greatly differ from a rectangular wing (ie: it will need to be high for extreme forward sweep angles). If the solution does not converge, increase this number as a first resort to try to force conversion
Minimum_Planform_Chord=0.00000000000000000001
Maximum_Planform_Chord=2
Mass=5
Gravitational_Constant=9.8055
Rho=1.1685
Viscosity=0.000017735

# PART 4 CONSTANTS
Washout_Adjust_Factor_Deg=0.00001 # Accurate to 0.001 degrees
J2=15000 # Note that this needs to be larger for span loadings that greatly differ from a rectangular wing (ie: it will need to be high for extreme forward sweep angles). If the solution does not converge, increase this number as a first resort to try to force conversion
Maximum_Negative_Twist_Deg=-20
Maximum_Positive_Twist_Deg=20

# PART 5 CONSTANTS
# Nil

# PART 6 CONSTANTS
Minimum_Structural_Chord=0.0977198697068404
Tip_Cap_Index=0.1 # Percentage of planform that is the tip
Smoothing_Addition_Root_Index1=0.2 # Percentage of planform where the smoothing starts at the tip (must be greater than zero)
Smoothing_Addition_Tip_Index1=2/pi # Percentage of planform where the smoothing ends at the root (must be less than one)
Smoothing_Addition_Root_Index2=2/pi # Percentage of planform where the smoothing starts at the tip (must be greater than zero)
Smoothing_Addition_Tip_Index2=0.99 # Percentage of planform where the smoothing ends at the root (must be less than one)
Fuselage_Length_Sizing_Factor=0.48493127382530354
Fuselage_Width_Sizing_Factor=1800
Fuselage_x_Offset_Index=0.1 # Percentage of fuselage being set back from the quarter chord line (ie: setting to 0.1 moves the quarter chord of the fuselage back by 10%)
Cm_Alpha_Root=-0.1178
Cm_Alpha_Tip=-0.001

#
#
#
#
#
#
#
#
#
#

# START OF CODE

#
#
#
#
#
#
#
#
#
#
# PART 1 - FUNCTION TO SOLVE THE WEISSINGER L-METHOD EQUATIONS GIVEN WING PLANFORM
# PART 1 - FUNCTION TO SOLVE THE WEISSINGER L-METHOD EQUATIONS GIVEN WING PLANFORM
# PART 1 - FUNCTION TO SOLVE THE WEISSINGER L-METHOD EQUATIONS GIVEN WING PLANFORM

# CONSTRAINTS ARE AS FOLLOWS
# 1. Only calculates symmetrical spanloads about the root chord (this is important as there are NaNs and Infs in the code that will present themselves in the LStar_v_Mu computation if an asymmetric loading is attempted)
# 2. No discontinuities in washout are permissable
# 3. The quarter chord line must be straight across the semispan
# 4. The wing must be planar
# 5. The flow is assumed to be inviscid (inapppropriate for flows which are dominated by viscous forces or when the the wing begins to stall)
# 6. There is only one chordwise control point (inappropriate for low aspect ratio wings (less than approximately AR=6) or for highly swept wings (sweep angles greater than 45 degrees))
# 7. There is an inbuilt error associated with approximation using Fourier transformations that scale with 1/(m^2)

function Weissinger_LMethod(FreeStream_Velocity,Sigma_SpanRatio,m,b,Root_Chord,Quarter_Chord_Sweep_Angle_Deg,Mach_Number,CL_RealWorld,Design_AoA_Deg,AlphaL0_Root_Deg,AlphaL0_Tip_Deg,Geometric_Twist_Distribution_Deg,Plotter_Step_Size)

# Preliminary calculations
m_Symm=Int64((m+1)/2)
Beta=sqrt(1-Mach_Number^2)
kv=CL_RealWorld/(2*pi/Beta)
Quarter_Chord_Sweep_Angle_Deg=atan(tan(Quarter_Chord_Sweep_Angle_Deg*pi/180)/Beta)*180/pi
Quarter_Chord_Sweep_Angle_Rad=Quarter_Chord_Sweep_Angle_Deg*pi/180

# Initial guess (obviously incorrect) chord distribution. Assume a rectangular wing of length b and chord uniformly equal to the root chord
c_v=ones(m_Symm,m_Symm).*Root_Chord
AR_v_1D=ones(m_Symm).*(b./c_v)*Beta*(1/kv)
AR_v_2D=ones(m_Symm,m_Symm).*AR_v_1D

# Setting up length vectors (default is m=7 for a 4 point symmetrical distribution)
Nu=1:1:m_Symm
v=1:1:m_Symm
n=1:1:m_Symm
Mu=0:1:m_Symm-1
Mu_1=1:2:m

# Making the above counter variables into 2D matrices for follow on calculations
v_Matrix_2D=(ones(m_Symm,m_Symm)).*v
n_Matrix_2D=(ones(m_Symm,m_Symm)).*n
n_Matrix_2DT=(ones(m_Symm,m_Symm)).*n.'
Mu_Matrix=ones(1,m_Symm).*(0:1:m_Symm-1) # Constructing Mu Matrix for fBar and LStar calculations

# Constructing Mu_1 Matrix in 3D for summation of f_n_Mu
Mu_1_Matrix=ones(1,1,m_Symm).*(0:2:m).+(ones(m_Symm,m_Symm,m_Symm)).*Mu_1.'
Mu_1_Matrix[findin(Mu_1_Matrix,Mu_1_Matrix[Mu_1_Matrix.>m])]=Mu_1_Matrix[findin(Mu_1_Matrix,Mu_1_Matrix[Mu_1_Matrix.>m])]-(m+1)
Mu_1_Matrix=permutedims(Mu_1_Matrix,[3,2,1])

# Polar coordinate transformations
Phi_v_2D=(ones(m_Symm,m_Symm)).*v_Matrix_2D.*pi/(m+1)
Phi_n_2D=(ones(m_Symm,m_Symm)).*n_Matrix_2D.*pi/(m+1)
Phi_n_2DT=(ones(m_Symm,m_Symm)).*n_Matrix_2DT.*pi/(m+1)
Phi_Mu=(ones(m_Symm,m_Symm)).*Mu_Matrix.*pi/(m+1)

# Nu position using polar coordinate tranformations
Nu_v=(ones(m_Symm,m_Symm)).*(cos.(Phi_v_2D))
NuBar_Mu=(ones(m_Symm,m_Symm)).*(cos.(Phi_Mu)).'

# Calculating a linear blend of airfoil from tip to root
Delta_AlphaL0_Deg=AlphaL0_Root_Deg-AlphaL0_Tip_Deg # Assumes a linear blend of airfoil across the wing in the following calculation
if Delta_AlphaL0_Deg==0
AlphaL0_Distribution_Deg=zeros(m_Symm,1)
else
AlphaL0_Distribution_Deg=(AlphaL0_Tip_Deg-AlphaL0_Root_Deg).*Nu_v[:,1]+(AlphaL0_Tip_Deg-(AlphaL0_Tip_Deg-AlphaL0_Root_Deg).*Nu_v[1,1])
end

# Calculation of angles
Design_AoA_Rad=(Design_AoA_Deg*pi/180).*ones(m_Symm,1)
AlphaL0_Distribution_Rad=(AlphaL0_Distribution_Deg*pi/180)
Geometric_Twist_Distribution_Rad=(Geometric_Twist_Distribution_Deg*pi/180).*ones(m_Symm,1)
Alpha_i=(1/(Sigma_SpanRatio^3))*(4.5*Sigma_SpanRatio-4-3*pi*(Sigma_SpanRatio-1).*Nu_v[:,1])
Alpha_i_Rad=Alpha_i/FreeStream_Velocity
Alpha=Design_AoA_Rad.-AlphaL0_Distribution_Rad.+Geometric_Twist_Distribution_Rad-Alpha_i_Rad # Note the plus on the geometric distribution. This is due to the method treating geometric distribution as an offset to the zero lift angle of attack (just as a negative sign on a cambered airfoil lowers the zero lift angle of attack). Therefore a negative geometric twist (normally would be washout) would increase lift which is counterintuitive without the plus sign on the geometric twist here

# Creation of additional variables to simplify following calculations
n_Additional_Matrix=((ones(m_Symm,m_Symm)).*(m+1-n)).'
n_Additional_Matrix[:,m_Symm]=0
Phi_n_Additional=(ones(m_Symm,m_Symm)).*n_Additional_Matrix.*pi/(m+1)

# Calculation of b_v_n, b_v_v and B_v_n
# Note that this is v rows x n columns
b_v_n=(sin.(Phi_n_2DT)./((cos.(Phi_n_2DT)-cos.(Phi_v_2D)).^2)).*((1-((-1).^(n_Matrix_2DT-v_Matrix_2D)))./(2*(m+1)))
b_v_v=eye(m_Symm,m_Symm).*(m+1)./(4*(sin.(Phi_v_2D)))
b_v_n[isnan(b_v_n)]=0 # Replacing NaN with zero
b_v_n[isinf(b_v_n)]=0 # Replacing Inf with zero
b_v_n_Additional=(sin.(Phi_n_Additional)./((cos.(Phi_n_Additional)-cos.(Phi_v_2D)).^2)).*((1-((-1).^(n_Additional_Matrix-v_Matrix_2D)))./(2*(m+1)))
B_v_n=b_v_n+b_v_n_Additional-b_v_v

# Calculation of f_n_Mu and fBar_n_Mu
# Note that this is n rows x Mu columns x Mu_1 pages (until permuted)
f_n_Mu=(2/(m+1)).*sum((Mu_1_Matrix.*sin.(Mu_1_Matrix.*Phi_n_2DT).*cos.(Mu_1_Matrix.*Phi_Mu)),3)
f_n_Mu=(f_n_Mu[:,:].').*ones(m_Symm,m_Symm,m_Symm)
fBar_n_Mu=permutedims(f_n_Mu,[3,1,2]).*2
fBar_n_Mu[:,:,1]=fBar_n_Mu[:,:,1]./2
fBar_n_Mu[:,m_Symm,:]=fBar_n_Mu[:,m_Symm,:]./2

# Calculation of LStar_v_Mu
# Note that this is v rows x Mu columns (until permuted)
K_Sub=Nu_v-NuBar_Mu
K_Add=Nu_v+NuBar_Mu
LStar_v_Mu_T1=(1./(AR_v_2D.*K_Sub)).*((sqrt(((1+AR_v_2D.*K_Sub.*tan(Quarter_Chord_Sweep_Angle_Rad)).^2)+((AR_v_2D).^2).*((K_Sub).^2)))-1)
LStar_v_Mu_T1[isnan(LStar_v_Mu_T1)]=tan(Quarter_Chord_Sweep_Angle_Rad) # Replacing NaN with tan(Quarter_Chord_Sweep_Angle_Rad) as this is the value the limit approached when Mu=v
LStar_v_Mu_T1[isinf(LStar_v_Mu_T1)]=tan(Quarter_Chord_Sweep_Angle_Rad) # Replacing Inf with tan(Quarter_Chord_Sweep_Angle_Rad) as this is the value the limit approached when Mu=v
LStar_v_Mu_T2=(1./(AR_v_2D.*K_Add)).*(((sqrt(((1+AR_v_2D.*K_Sub.*tan(Quarter_Chord_Sweep_Angle_Rad)).^2)+((AR_v_2D).^2).*((K_Add).^2)))./(1+2.*AR_v_2D.*Nu_v.*tan(Quarter_Chord_Sweep_Angle_Rad)))-1)
LStar_v_Mu_T3=2.*tan(Quarter_Chord_Sweep_Angle_Rad).*(sqrt(((1+AR_v_2D.*Nu_v.*tan(Quarter_Chord_Sweep_Angle_Rad)).^2)+((AR_v_2D.^2).*(Nu_v.^2)))./(1+2*AR_v_2D.*Nu_v.*tan(Quarter_Chord_Sweep_Angle_Rad)))
LStar_v_Mu=(LStar_v_Mu_T1-LStar_v_Mu_T2-LStar_v_Mu_T3).*ones(m_Symm,m_Symm,m_Symm)
LStar_v_Mu=permutedims(LStar_v_Mu,[1,3,2])

# Calculation of gBar_v_n
gBar_v_n=(-1/(2*(m+1))).*sum(fBar_n_Mu.*LStar_v_Mu,3)
gBar_v_n=gBar_v_n[:,:]

# Calculation of a_v_n
a_v_n=-2*B_v_n.+(AR_v_1D.*gBar_v_n)

# Solve for G_n_a
Gna=a_v_n\Alpha

# Re-dimensionalize
Gna=Gna*b*FreeStream_Velocity
#Gna=Gna/Gna[m_Symm]

# PLOTTER
Nu_Line=0:Plotter_Step_Size:1
Total_Plot_Points=length(Nu_Line)
v_Matrix_2D_mPoint=(ones(m,m)).*(1:1:m)
G_Phi=zeros(Total_Plot_Points,1)
Gna_Plotter=cat(1,Gna[1:length(Gna)-1],flipdim(Gna,1))
for i=1:Total_Plot_Points
Nu_Plotter=zeros(m,m)+Plotter_Step_Size*i-Plotter_Step_Size
Phi_Plotter=acos(Nu_Plotter)
Phi_Plotter_Matrix=sin.(v_Matrix_2D_mPoint.*((v_Matrix_2D_mPoint.*pi/(m+1)).')).*sin(v_Matrix_2D_mPoint.*Phi_Plotter)
G_Phi_Matrix=sum(Phi_Plotter_Matrix*Gna_Plotter,1)
G_Phi[i,1]=(2/(m+1)).*G_Phi_Matrix[1]
end

Nu_v_OUTPUT=Nu_v[:,1]
c_v_OUTPUT=c_v[:,1]
return Nu_v_OUTPUT,c_v_OUTPUT,Nu_Line,G_Phi,Gna,Alpha,Geometric_Twist_Distribution_Deg
end

# END PART 1 - FUNCTION TO SOLVE THE WEISSINGER L-METHOD EQUATIONS GIVEN WING PLANFORM
# END PART 1 - FUNCTION TO SOLVE THE WEISSINGER L-METHOD EQUATIONS GIVEN WING PLANFORM
# END PART 1 - FUNCTION TO SOLVE THE WEISSINGER L-METHOD EQUATIONS GIVEN WING PLANFORM
#
#
#
#
#
#
#
#
#
#
# PART 2 - FUNCTION TO GENERATE m VALUE FOR GRID INDEPENDENCE
# PART 2 - FUNCTION TO GENERATE m VALUE FOR GRID INDEPENDENCE
# PART 2 - FUNCTION TO GENERATE m VALUE FOR GRID INDEPENDENCE

function Grid_Independence_Iterator(FreeStream_Velocity,Sigma_SpanRatio,m,b,Root_Chord,Quarter_Chord_Sweep_Angle_Deg,Mach_Number,CL_RealWorld,Design_AoA_Deg,AlphaL0_Root_Deg,AlphaL0_Tip_Deg,Geometric_Twist_Distribution_Deg,Plotter_Step_Size,m_Limit_Computer,m_Increment,Accuracy_Increment)

# Using m=599 (300 point symmetrical calculation) due to computational storage limitations
Nu_v_GRID,c_v_GRID,Nu_Line_GRID,G_Phi_GRID,Gna_GRID,Alpha_GRID,Geometric_Twist_Distribution_Deg_GRID=Weissinger_LMethod(FreeStream_Velocity,Sigma_SpanRatio,m_Limit_Computer,b,Root_Chord,Quarter_Chord_Sweep_Angle_Deg,Mach_Number,CL_RealWorld,Design_AoA_Deg,AlphaL0_Root_Deg,AlphaL0_Tip_Deg,Geometric_Twist_Distribution_Deg,Plotter_Step_Size)

# Initialisation
println(" ")
println("COMMENCING GRID INDEPENDENCE STUDY")
println(" ")
Accuracy_Metric=1
m_ITER=m-m_Increment
Nu_v_ACTUAL,c_v_ACTUAL,Nu_Line_ACTUAL,G_Phi_ACTUAL,Gna_ACTUAL,Alpha_ACTUAL,Geometric_Twist_Distribution_Deg_ACTUAL=Weissinger_LMethod(FreeStream_Velocity,Sigma_SpanRatio,m_ITER,b,Root_Chord,Quarter_Chord_Sweep_Angle_Deg,Mach_Number,CL_RealWorld,Design_AoA_Deg,AlphaL0_Root_Deg,AlphaL0_Tip_Deg,Geometric_Twist_Distribution_Deg,Plotter_Step_Size)

Accuracy_Metric_Singular=1
while Accuracy_Metric_Singular>Accuracy_Increment
m_ITER=m_ITER+m_Increment

# Using the actual m inputted
Nu_v_ACTUAL,c_v_ACTUAL,Nu_Line_ACTUAL,G_Phi_ACTUAL,Gna_ACTUAL,Alpha_ACTUAL,Geometric_Twist_Distribution_Deg_ACTUAL=Weissinger_LMethod(FreeStream_Velocity,Sigma_SpanRatio,m_ITER,b,Root_Chord,Quarter_Chord_Sweep_Angle_Deg,Mach_Number,CL_RealWorld,Design_AoA_Deg,AlphaL0_Root_Deg,AlphaL0_Tip_Deg,Geometric_Twist_Distribution_Deg,Plotter_Step_Size)

Accuracy_Metric=abs(G_Phi_GRID.-G_Phi_ACTUAL)
Accuracy_Metric_Singular=maximum(Accuracy_Metric)

println("Percantage Error= ",Accuracy_Metric_Singular*100," %")
println("m = ",m_ITER)
end

return Nu_v_ACTUAL,c_v_ACTUAL,m_ITER,Gna_ACTUAL,Nu_Line_GRID,G_Phi_GRID,Nu_Line_ACTUAL,G_Phi_ACTUAL,Alpha_ACTUAL,Geometric_Twist_Distribution_Deg_ACTUAL,Accuracy_Metric
end

# END PART 2 - FUNCTION TO GENERATE m VALUE FOR GRID INDEPENDENCE
# END PART 2 - FUNCTION TO GENERATE m VALUE FOR GRID INDEPENDENCE
# END PART 2 - FUNCTION TO GENERATE m VALUE FOR GRID INDEPENDENCE
#
#
#
#
#
#
#
#
#
#
# PART 3 - FUNCTION TO CALCULATE REQUIRED PLANFORM
# PART 3 - FUNCTION TO CALCULATE REQUIRED PLANFORM
# PART 3 - FUNCTION TO CALCULATE REQUIRED PLANFORM

function Planform_Optimizer(FreeStream_Velocity,Sigma_SpanRatio,m,b,Root_Chord,Quarter_Chord_Sweep_Angle_Deg,Mach_Number,CL_RealWorld,Design_AoA_Deg,AlphaL0_Root_Deg,AlphaL0_Tip_Deg,Geometric_Twist_Distribution_Deg,Plotter_Step_Size,m_Limit_Computer,m_Increment,Accuracy_Increment,Planform_Adjust_Factor,J1,Minimum_Planform_Chord,Maximum_Planform_Chord,Mass,Gravitational_Constant,Rho,Viscosity)

println(" ")
println("COMMENCING PLANFORM OPTIMIZATION")
println(" ")

# Initializing using a zeroed start state planform
m_Symm_IND=Int64((m_IND+1)/2)
c_v_STARTSTATE=ones(m_Symm_IND,1).*Minimum_Planform_Chord # Getting rid of tip errors by starting from the minimum planform
c_v_STARTSTATE[m_Symm_IND,1]=Root_Chord # Getting rid of tip errors by starting from the minimum planform
Nu_v_PLANFORM,c_v_PLANFORM,Nu_Line_PLANFORM,G_Phi_PLANFORM,Gna_PLANFORM,Alpha_PLANFORM,Geometric_Twist_Distribution_Deg_PLANFORM=Weissinger_LMethod(FreeStream_Velocity,Sigma_SpanRatio,m_IND,b,c_v_STARTSTATE,Quarter_Chord_Sweep_Angle_Deg,Mach_Number,CL_RealWorld,Design_AoA_Deg,AlphaL0_Root_Deg,AlphaL0_Tip_Deg,Geometric_Twist_Distribution_Deg,Plotter_Step_Size)

# Preliminary calculations
Nu_v_PLANFORM=abs(Nu_v_PLANFORM) #Avoiding slighly negative number near zero which becomes problematic when computing asech function
Gamma_BSLD=sqrt(1-Nu_v_PLANFORM.^2)+6*((1-Sigma_SpanRatio)/(3*Sigma_SpanRatio-2))*(Nu_v_PLANFORM.^2).*asech(Nu_v_PLANFORM)
Gamma_BSLD=Gamma_BSLD*Mass*Gravitational_Constant/(2*b*Rho*FreeStream_Velocity)
Planform_Upper_Limit=Gamma_BSLD+ones(m_Symm_IND,1).*Planform_Adjust_Factor
Planform_Lower_Limit=Gamma_BSLD-ones(m_Symm_IND,1).*Planform_Adjust_Factor

# Transforming the dimensional circulation from the Weissinger function into Lift (per unit span)
#Gna_PLANFORM=Gna_PLANFORM*Rho*FreeStream_Velocity

# Start of loop to iterate to a converged solution
j=0
for k=1:m_Symm_IND
while (Gna_PLANFORM[k,1]<Planform_Lower_Limit[k,1]||Gna_PLANFORM[k,1]>Planform_Upper_Limit[k,1])&&j<J1
j=j+1
for i=1:m_Symm_IND
if Gna_PLANFORM[i,1]<Gamma_BSLD[i,1]
c_v_PLANFORM_INCREMENT=(abs(Gamma_BSLD[i,1])-abs(Gna_PLANFORM[i,1]))*abs(Gamma_BSLD[i,1])*cos(Quarter_Chord_Sweep_Angle_Deg*pi/180) # The converging increment scales with Gamma value (large at root and small at tip) as well as the cosine of the sweep angle (large sweep angles converge more slowly to ensure the solution converges)
c_v_PLANFORM[i,1]=c_v_PLANFORM[i,1]+c_v_PLANFORM_INCREMENT
if c_v_PLANFORM[i,1]>Maximum_Planform_Chord
c_v_PLANFORM[i,1]=Maximum_Planform_Chord
end
end
if Gna_PLANFORM[i,1]>Gamma_BSLD[i,1]
c_v_PLANFORM_INCREMENT=(abs(Gna_PLANFORM[i,1])-abs(Gamma_BSLD[i,1]))*abs(Gamma_BSLD[i,1])*cos(Quarter_Chord_Sweep_Angle_Deg*pi/180) # The converging increment scales with Gamma value (large at root and small at tip) as well as the cosine of the sweep angle (large sweep angles converge more slowly to ensure the solution converges)
c_v_PLANFORM[i,1]=c_v_PLANFORM[i,1]-c_v_PLANFORM_INCREMENT
if c_v_PLANFORM[i,1]<Minimum_Planform_Chord
c_v_PLANFORM[i,1]=Minimum_Planform_Chord
end
end
end

# Calling Weissinger function to optimize planform
Nu_v_PLANFORM,c_v_PLANFORM,Nu_Line_PLANFORM,G_Phi_PLANFORM,Gna_PLANFORM,Alpha_PLANFORM,Geometric_Twist_Distribution_Deg_PLANFORM=Weissinger_LMethod(FreeStream_Velocity,Sigma_SpanRatio,m_IND,b,c_v_PLANFORM,Quarter_Chord_Sweep_Angle_Deg,Mach_Number,CL_RealWorld,Design_AoA_Deg,AlphaL0_Root_Deg,AlphaL0_Tip_Deg,Geometric_Twist_Distribution_Deg,Plotter_Step_Size)

# Displaying residuals (1.0 for convergence criteria met)
Percentage_Residuals=Planform_Adjust_Factor./abs(Gna_PLANFORM[:,1]-Gamma_BSLD[:,1])
Percentage_Residuals[findin(Percentage_Residuals,Percentage_Residuals[Percentage_Residuals.>1])]=1
Percentage_Residuals[isinf(Percentage_Residuals)]=1
Percentage_Residuals=1-Percentage_Residuals
println("j1 = ",j," out of ",J1)
println(" ")
println("Residuals = ",Percentage_Residuals)
println(" ")
println("Highest Residual = ",maximum(Percentage_Residuals))
println("Average Residual = ",mean(Percentage_Residuals))
println(" ")
end
end

println(" ")
println("PLANFORM OPTIMIZATION COMPLETE")
println(" ")

# Final residuals
Nu_v_PLANFORM,c_v_PLANFORM,Nu_Line_PLANFORM,G_Phi_PLANFORM,Gna_PLANFORM,Alpha_PLANFORM,Geometric_Twist_Distribution_Deg_PLANFORM=Weissinger_LMethod(FreeStream_Velocity,Sigma_SpanRatio,m_IND,b,c_v_PLANFORM,Quarter_Chord_Sweep_Angle_Deg,Mach_Number,CL_RealWorld,Design_AoA_Deg,AlphaL0_Root_Deg,AlphaL0_Tip_Deg,Geometric_Twist_Distribution_Deg,Plotter_Step_Size)
Percentage_Residuals=Planform_Adjust_Factor./abs(Gna_PLANFORM[:,1]-Gamma_BSLD[:,1])
Percentage_Residuals[findin(Percentage_Residuals,Percentage_Residuals[Percentage_Residuals.>1])]=1
Percentage_Residuals[isinf(Percentage_Residuals)]=1
Percentage_Residuals=1-Percentage_Residuals
println("Final Residuals = ",Percentage_Residuals)
println(" ")
println("Final Highest Residual = ",maximum(Percentage_Residuals))
println(" ")
println("Alpha (Degrees) = ",Alpha_PLANFORM*180/pi)
println(" ")

# Scaled plotting of the planform
Nu_v_PLANFORM_SCALED=Nu_v_PLANFORM*b/2
Quarter_Chord_Sweep_Line=Nu_v_PLANFORM_SCALED.*tan(-Quarter_Chord_Sweep_Angle_Deg*pi/180)
LE_Line=Quarter_Chord_Sweep_Line+0.25*c_v_PLANFORM
TE_Line=Quarter_Chord_Sweep_Line-0.75*c_v_PLANFORM

# Find the AR of the final wing
Total_Area=0
for i=1:m_Symm_IND-1
y_Distance=abs(Nu_v_PLANFORM_SCALED[i+1]-Nu_v_PLANFORM_SCALED[i])
x_Distance=(c_v_PLANFORM[i+1]+c_v_PLANFORM[i])/2
Area=y_Distance*x_Distance
Total_Area=Total_Area+Area
end

Aspect_Ratio=(b^2)/(Total_Area*2)
println("Aspect Ratio = ",Aspect_Ratio)

Reynolds_Number_Root=Rho*FreeStream_Velocity*c_v_PLANFORM[m_Symm_IND]/Viscosity
println("Root Reynolds Number = ",Reynolds_Number_Root/1000," K ")

Reynolds_Number_PseudoAverage=mean(Rho*FreeStream_Velocity*c_v_PLANFORM/Viscosity)
println("Pseudo-Average Reynolds Number = ",Reynolds_Number_PseudoAverage/1000," K ")

return Nu_Line_GRID,Nu_Line_ACTUAL,Nu_Line_PLANFORM,G_Phi_GRID,G_Phi_ACTUAL,Accuracy_Metric,Nu_v_PLANFORM,c_v_PLANFORM,Gamma_BSLD,G_Phi_PLANFORM,Nu_v_PLANFORM_SCALED,Quarter_Chord_Sweep_Line,LE_Line,TE_Line,Aspect_Ratio,m_Symm_IND,Alpha_PLANFORM,Geometric_Twist_Distribution_Deg_PLANFORM
end

# END PART 3 - FUNCTION TO CALCULATE REQUIRED PLANFORM
# END PART 3 - FUNCTION TO CALCULATE REQUIRED PLANFORM
# END PART 3 - FUNCTION TO CALCULATE REQUIRED PLANFORM
#
#
#
#
#
#
#
#
#
#
# PART 4 - FUNCTION TO CALCULATE REQUIRED WASHOUT
# PART 4 - FUNCTION TO CALCULATE REQUIRED WASHOUT
# PART 4 - FUNCTION TO CALCULATE REQUIRED WASHOUT

function Washout_Optimizer(FreeStream_Velocity,Sigma_SpanRatio,m,b,Root_Chord,Quarter_Chord_Sweep_Angle_Deg,Mach_Number,CL_RealWorld,Design_AoA_Deg,AlphaL0_Root_Deg,AlphaL0_Tip_Deg,Geometric_Twist_Distribution_Deg,Plotter_Step_Size,m_Limit_Computer,m_Increment,Accuracy_Increment,Washout_Adjust_Factor_Deg,J2,Minimum_Planform_Chord,Maximum_Planform_Chord,Mass,Gravitational_Constant,Rho,Viscosity)

println(" ")
println("COMMENCING WASHOUT OPTIMIZATION")
println(" ")

# Initializing using a rectangular planform
m_Symm_IND=Int64((m_IND+1)/2)
Geometric_Twist_Distribution_Deg_STARTSRATE=zeros(m_Symm_IND,1)
Nu_v_WASHOUT,c_v_WASHOUT,Nu_Line_WASHOUT,G_Phi_WASHOUT,Gna_WASHOUT,Alpha_WASHOUT,Geometric_Twist_Distribution_Deg_WASHOUT=Weissinger_LMethod(FreeStream_Velocity,Sigma_SpanRatio,m_IND,b,Root_Chord,Quarter_Chord_Sweep_Angle_Deg,Mach_Number,CL_RealWorld,Design_AoA_Deg,AlphaL0_Root_Deg,AlphaL0_Tip_Deg,Geometric_Twist_Distribution_Deg_STARTSRATE,Plotter_Step_Size)

# Preliminary calculations
Washout_Adjust_Factor_Rad=Washout_Adjust_Factor_Deg*pi/180
Nu_v_WASHOUT=abs(Nu_v_WASHOUT) #Avoiding slighly negative number near zero which becomes problematic when computing asech function
Gamma_BSLD=sqrt(1-Nu_v_WASHOUT.^2)+6*((1-Sigma_SpanRatio)/(3*Sigma_SpanRatio-2))*(Nu_v_WASHOUT.^2).*asech(Nu_v_WASHOUT)
Gamma_BSLD=Gamma_BSLD*Mass*Gravitational_Constant/(2*b*Rho*FreeStream_Velocity)
Washout_Upper_Limit=Gamma_BSLD+ones(m_Symm_IND,1).*Washout_Adjust_Factor_Rad
Washout_Lower_Limit=Gamma_BSLD-ones(m_Symm_IND,1).*Washout_Adjust_Factor_Rad

# Transforming the dimensional circulation from the Weissinger function into Lift (per unit span)
#Gna_WASHOUT=Gna_WASHOUT*Rho*FreeStream_Velocity

# Start of loop to iterate to a converged solution
j=0
for k=1:m_Symm_IND
while (Gna_WASHOUT[k,1]<Washout_Lower_Limit[k,1]||Gna_WASHOUT[k,1]>Washout_Upper_Limit[k,1])&&j<J2
j=j+1
for i=1:m_Symm_IND
if Gna_WASHOUT[i,1]<Gamma_BSLD[i,1]
Geometric_Twist_Distribution_Deg_INCREMENT=(Gamma_BSLD[i,1]-Gna_WASHOUT[i,1])*cos(Quarter_Chord_Sweep_Angle_Deg*pi/180)*(Nu_v_WASHOUT[i,1]+1) # The converging increment scales with the cosine of the sweep angle (large sweep angles converge more slowly to ensure the solution converges)
Geometric_Twist_Distribution_Deg_WASHOUT[i,1]=Geometric_Twist_Distribution_Deg_WASHOUT[i,1]+Geometric_Twist_Distribution_Deg_INCREMENT
if Geometric_Twist_Distribution_Deg_WASHOUT[i,1]>Maximum_Positive_Twist_Deg
Geometric_Twist_Distribution_Deg_WASHOUT[i,1]=Maximum_Positive_Twist_Deg
end
end
if Gna_WASHOUT[i,1]>Gamma_BSLD[i,1]
Geometric_Twist_Distribution_Deg_INCREMENT=(Gna_WASHOUT[i,1]-Gamma_BSLD[i,1])*cos(Quarter_Chord_Sweep_Angle_Deg*pi/180)*(Nu_v_WASHOUT[i,1]+1) # The converging increment scales with the cosine of the sweep angle (large sweep angles converge more slowly to ensure the solution converges)
Geometric_Twist_Distribution_Deg_WASHOUT[i,1]=Geometric_Twist_Distribution_Deg_WASHOUT[i,1]-Geometric_Twist_Distribution_Deg_INCREMENT
if Geometric_Twist_Distribution_Deg_WASHOUT[i,1]<Maximum_Negative_Twist_Deg
Geometric_Twist_Distribution_Deg_WASHOUT[i,1]=Maximum_Negative_Twist_Deg
end
end
end

# Calling Weissinger function to optimize planform
Nu_v_WASHOUT,c_v_WASHOUT,Nu_Line_WASHOUT,G_Phi_WASHOUT,Gna_WASHOUT,Alpha_WASHOUT,Geometric_Twist_Distribution_Deg_WASHOUT=Weissinger_LMethod(FreeStream_Velocity,Sigma_SpanRatio,m_IND,b,Root_Chord,Quarter_Chord_Sweep_Angle_Deg,Mach_Number,CL_RealWorld,Design_AoA_Deg,AlphaL0_Root_Deg,AlphaL0_Tip_Deg,Geometric_Twist_Distribution_Deg_WASHOUT,Plotter_Step_Size)

# Displaying residuals (1.0 for convergence criteria met)
Percentage_Residuals=Washout_Adjust_Factor_Rad./abs(Gna_WASHOUT[:,1]-Gamma_BSLD[:,1])
Percentage_Residuals[findin(Percentage_Residuals,Percentage_Residuals[Percentage_Residuals.>1])]=1
Percentage_Residuals[isinf(Percentage_Residuals)]=1
Percentage_Residuals=1-Percentage_Residuals
println("j2 = ",j," out of ",J2)
println(" ")
println("Residuals = ",Percentage_Residuals)
println(" ")
println("Highest Residual = ",maximum(Percentage_Residuals))
println("Average Residual = ",mean(Percentage_Residuals))
println(" ")
end
end

println(" ")
println("WASHOUT OPTIMIZATION COMPLETE")
println(" ")

# Final residuals
Nu_v_WASHOUT,c_v_WASHOUT,Nu_Line_WASHOUT,G_Phi_WASHOUT,Gna_WASHOUT,Alpha_WASHOUT,Geometric_Twist_Distribution_Deg_WASHOUT=Weissinger_LMethod(FreeStream_Velocity,Sigma_SpanRatio,m_IND,b,Root_Chord,Quarter_Chord_Sweep_Angle_Deg,Mach_Number,CL_RealWorld,Design_AoA_Deg,AlphaL0_Root_Deg,AlphaL0_Tip_Deg,Geometric_Twist_Distribution_Deg_WASHOUT,Plotter_Step_Size)
Percentage_Residuals=Washout_Adjust_Factor_Rad./abs(Gna_WASHOUT[:,1]-Gamma_BSLD[:,1])
Percentage_Residuals[findin(Percentage_Residuals,Percentage_Residuals[Percentage_Residuals.>1])]=1
Percentage_Residuals[isinf(Percentage_Residuals)]=1
Percentage_Residuals=1-Percentage_Residuals
println("Final Residuals = ",Percentage_Residuals)
println(" ")
println("Final Highest Residual = ",maximum(Percentage_Residuals))
println(" ")
println("Alpha (Degrees) = ",Alpha_WASHOUT*180/pi)
println(" ")

# Scaled plotting of the planform
Nu_v_WASHOUT_SCALED=Nu_v_WASHOUT*b/2
Quarter_Chord_Sweep_Line=Nu_v_WASHOUT_SCALED.*tan(-Quarter_Chord_Sweep_Angle_Deg*pi/180)
LE_Line=Quarter_Chord_Sweep_Line+0.25*c_v_WASHOUT
TE_Line=Quarter_Chord_Sweep_Line-0.75*c_v_WASHOUT

# Find the AR of the final wing
Total_Area=0
for i=1:m_Symm_IND-1
y_Distance=abs(Nu_v_WASHOUT_SCALED[i+1]-Nu_v_WASHOUT_SCALED[i])
x_Distance=(c_v_WASHOUT[i+1]+c_v_WASHOUT[i])/2
Area=y_Distance*x_Distance
Total_Area=Total_Area+Area
end

Aspect_Ratio=(b^2)/(Total_Area*2)
println("Aspect Ratio = ",Aspect_Ratio)

Reynolds_Number_Root=Rho*FreeStream_Velocity*c_v_WASHOUT[m_Symm_IND]/Viscosity
println("Root Reynolds Number = ",Reynolds_Number_Root/1000," K ")

Reynolds_Number_PseudoAverage=mean(Rho*FreeStream_Velocity*c_v_WASHOUT/Viscosity)
println("Pseudo-Average Reynolds Number = ",Reynolds_Number_PseudoAverage/1000," K ")

return Nu_Line_GRID,Nu_Line_ACTUAL,Nu_Line_WASHOUT,G_Phi_GRID,G_Phi_ACTUAL,Accuracy_Metric,Nu_v_WASHOUT,c_v_WASHOUT,Gamma_BSLD,G_Phi_WASHOUT,Nu_v_WASHOUT_SCALED,Quarter_Chord_Sweep_Line,LE_Line,TE_Line,Aspect_Ratio,m_Symm_IND,Alpha_WASHOUT,Geometric_Twist_Distribution_Deg_WASHOUT
end

# END PART 4 - FUNCTION TO CALCULATE REQUIRED WASHOUT
# END PART 4 - FUNCTION TO CALCULATE REQUIRED WASHOUT
# END PART 4 - FUNCTION TO CALCULATE REQUIRED WASHOUT
#
#
#
#
#
#
#
#
#
#
# PART 5 - INITIALISATION OF FINAL OPTIMIZATION
# PART 5 - INITIALISATION OF FINAL OPTIMIZATION
# PART 5 - INITIALISATION OF FINAL OPTIMIZATION

# Calling grid independence function
Nu_v_IND,c_v_IND,m_IND,Gna_IND,Nu_Line_GRID,G_Phi_GRID,Nu_Line_ACTUAL,G_Phi_ACTUAL,Alpha_ACTUAL,Geometric_Twist_Distribution_Deg_ACTUAL,Accuracy_Metric=Grid_Independence_Iterator(FreeStream_Velocity,Sigma_SpanRatio,m,b,Root_Chord,Quarter_Chord_Sweep_Angle_Deg,Mach_Number,CL_RealWorld,Design_AoA_Deg,AlphaL0_Root_Deg,AlphaL0_Tip_Deg,Geometric_Twist_Distribution_Deg,Plotter_Step_Size,m_Limit_Computer,m_Increment,Accuracy_Increment)

println(" ")
println("GRID INDEPENDENCE STUDY COMPLETE")
println("m = ",m_IND," for Grid Independence")
println(" ")

# CALLING THE PART 1, PART 2, PART 3 - PLANFORM OPTIMIZER
Nu_Line_GRID,Nu_Line_ACTUAL,Nu_Line_PLANFORM,G_Phi_GRID,G_Phi_ACTUAL,Accuracy_Metric,Nu_v_PLANFORM,c_v_PLANFORM,Gamma_BSLD,G_Phi_PLANFORM,Nu_v_PLANFORM_SCALED,Quarter_Chord_Sweep_Line,LE_Line,TE_Line,Aspect_Ratio,m_Symm_IND,Alpha_PLANFORM,Geometric_Twist_Distribution_Deg_PLANFORM=Planform_Optimizer(FreeStream_Velocity,Sigma_SpanRatio,m_IND,b,Root_Chord,Quarter_Chord_Sweep_Angle_Deg,Mach_Number,CL_RealWorld,Design_AoA_Deg,AlphaL0_Root_Deg,AlphaL0_Tip_Deg,Geometric_Twist_Distribution_Deg,Plotter_Step_Size,m_Limit_Computer,m_Increment,Accuracy_Increment,Planform_Adjust_Factor,J1,Minimum_Planform_Chord,Maximum_Planform_Chord,Mass,Gravitational_Constant,Rho,Viscosity)

# De-comment for plots of the optimal planform to achieve the Gamma Distribution without any washout required
#=
using Plots
plotly()
plot()
plot!(Nu_v_PLANFORM,c_v_PLANFORM)
plot!(Nu_v_PLANFORM,Gamma_BSLD)
plot!(Nu_Line_PLANFORM,G_Phi_PLANFORM)
gui()

using Plots
plotly()
plot()
xlims!(-1.1*b/2,1.1*b/2)
ylims!(-0.75*1.1*b/2,0.75*1.1*b/2)
plot!(Nu_v_PLANFORM_SCALED,Quarter_Chord_Sweep_Line)
plot!(Nu_v_PLANFORM_SCALED,LE_Line)
plot!(Nu_v_PLANFORM_SCALED,TE_Line)
plot!(-Nu_v_PLANFORM_SCALED,Quarter_Chord_Sweep_Line)
plot!(-Nu_v_PLANFORM_SCALED,LE_Line)
plot!(-Nu_v_PLANFORM_SCALED,TE_Line)
plot!(Nu_v_PLANFORM_SCALED,c_v_PLANFORM+1)
plot!(Nu_v_PLANFORM_SCALED,0.15.*ones(m_Symm_IND,1)+1) # 150mm mark for acceptable Reynolds Number
plot!(-Nu_v_PLANFORM_SCALED,c_v_PLANFORM+1)
plot!(-Nu_v_PLANFORM_SCALED,0.15.*ones(m_Symm_IND,1)+1) # 150mm mark for acceptable Reynolds Number
gui()

using Plots
plotly()
plot()
plot!(Nu_v_PLANFORM,Alpha_PLANFORM*180/pi)
plot!(Nu_v_PLANFORM,Geometric_Twist_Distribution_Deg_PLANFORM)
gui()
=#

# END PART 5 - INITIALISATION OF FINAL OPTIMIZATION
# END PART 5 - INITIALISATION OF FINAL OPTIMIZATION
# END PART 5 - INITIALISATION OF FINAL OPTIMIZATION
#
#
#
#
#
#
#
#
#
#
# PART 6 - FINAL OPTIMIZATION OF MODIFIED PLANFORM
# PART 6 - FINAL OPTIMIZATION OF MODIFIED PLANFORM
# PART 6 - FINAL OPTIMIZATION OF MODIFIED PLANFORM

m_Symm_IND=Int64((m_IND+1)/2)
c_v_PLANFORM_ADDITIONAL=zeros(m_Symm_IND)

# PLANFORM RULE NUMBER 1
# If the planform is less than the minimum strutural chord within the non tip cap component of the wing, then planform must be added (the tip cap index of the wing indicates the tip cap and this rule does not apply)
Planform_Builder=0
for l=1:m_Symm_IND
if c_v_PLANFORM[l]<Minimum_Structural_Chord&&Nu_v_PLANFORM[l]<=(1-Tip_Cap_Index)
c_v_PLANFORM_ADDITIONAL[l]=Minimum_Structural_Chord-c_v_PLANFORM[l]
Planform_Builder=1
end
end

# PLANFORM RULE NUMBER 2
# The planform of the tip cap of the wing must be convex (elliptical curve used)
Tip_Cap_StartLocation=findin(c_v_PLANFORM_ADDITIONAL,maximum(c_v_PLANFORM_ADDITIONAL))
Elliptical_x_Offset=Nu_v_PLANFORM[Tip_Cap_StartLocation]
Elliptical_x_Offset=Elliptical_x_Offset[1]
Elliptical_y_Offset=0
Elliptical_x_Radius=1-Elliptical_x_Offset
Elliptical_y_Offset=c_v_PLANFORM_ADDITIONAL[Tip_Cap_StartLocation]+c_v_PLANFORM[Tip_Cap_StartLocation]
Elliptical_y_Offset=Elliptical_y_Offset[1]
for l=1:m_Symm_IND
if Planform_Builder==1&&Nu_v_PLANFORM[l]>=(1-Tip_Cap_Index)
c_v_PLANFORM_ADDITIONAL[l]=sqrt((1-(((Nu_v_PLANFORM[l]-Elliptical_x_Offset)^2)/(Elliptical_x_Radius^2)))*((Elliptical_y_Offset)^2))
c_v_PLANFORM_ADDITIONAL[l]=c_v_PLANFORM_ADDITIONAL[l]-c_v_PLANFORM[l]
end
end

# PLANFORM RULE NUMBER 3
# The planform must have no sharp angles or discontinuities due to the addition of the minimum structural chord
Smoothing_Addition_Root_Index=Smoothing_Addition_Root_Index1
Smoothing_Addition_Tip_Index=Smoothing_Addition_Tip_Index1
Nu_v_PLANFORM_SMOOTHING=zeros(Nu_v_PLANFORM)+Nu_v_PLANFORM
Nu_v_PLANFORM_SMOOTHING[findin(Nu_v_PLANFORM_SMOOTHING,Nu_v_PLANFORM_SMOOTHING[Nu_v_PLANFORM_SMOOTHING.<=Smoothing_Addition_Root_Index])]=0
Nu_v_PLANFORM_SMOOTHING[findin(Nu_v_PLANFORM_SMOOTHING,Nu_v_PLANFORM_SMOOTHING[Nu_v_PLANFORM_SMOOTHING.>=Smoothing_Addition_Tip_Index])]=0
Smoothing_Addition_Root_IndexLocator=findfirst(Nu_v_PLANFORM_SMOOTHING)
Smoothing_Addition_Tip_IndexLocator=findlast(Nu_v_PLANFORM_SMOOTHING)
Smoothing_Addition_Tip_Gradient=((c_v_PLANFORM[Smoothing_Addition_Root_IndexLocator]+c_v_PLANFORM_ADDITIONAL[Smoothing_Addition_Root_IndexLocator])-(c_v_PLANFORM[Smoothing_Addition_Root_IndexLocator-1]+c_v_PLANFORM_ADDITIONAL[Smoothing_Addition_Root_IndexLocator-1]))/(Nu_v_PLANFORM[Smoothing_Addition_Root_IndexLocator]-Nu_v_PLANFORM[Smoothing_Addition_Root_IndexLocator-1])
Smoothing_Addition_Root_Gradient=((c_v_PLANFORM[Smoothing_Addition_Tip_IndexLocator+1]+c_v_PLANFORM_ADDITIONAL[Smoothing_Addition_Tip_IndexLocator+1])-(c_v_PLANFORM[Smoothing_Addition_Tip_IndexLocator]+c_v_PLANFORM_ADDITIONAL[Smoothing_Addition_Tip_IndexLocator]))/(Nu_v_PLANFORM[Smoothing_Addition_Tip_IndexLocator+1]-Nu_v_PLANFORM[Smoothing_Addition_Tip_IndexLocator])
Cubic_x_Tip=Nu_v_PLANFORM[Smoothing_Addition_Root_IndexLocator]
Cubic_x_Root=Nu_v_PLANFORM[Smoothing_Addition_Tip_IndexLocator]
Cubic_y_Tip=c_v_PLANFORM[Smoothing_Addition_Root_IndexLocator]+c_v_PLANFORM_ADDITIONAL[Smoothing_Addition_Root_IndexLocator]
Cubic_y_Root=c_v_PLANFORM[Smoothing_Addition_Tip_IndexLocator]+c_v_PLANFORM_ADDITIONAL[Smoothing_Addition_Tip_IndexLocator]

Cubic_Matrix=zeros(4,4)
Cubic_Matrix[1,1]=Cubic_x_Tip^3
Cubic_Matrix[2,1]=Cubic_x_Root^3
Cubic_Matrix[3,1]=3*Cubic_x_Tip^2
Cubic_Matrix[4,1]=3*Cubic_x_Root^2
Cubic_Matrix[1,2]=Cubic_x_Tip^2
Cubic_Matrix[2,2]=Cubic_x_Root^2
Cubic_Matrix[3,2]=2*Cubic_x_Tip
Cubic_Matrix[4,2]=2*Cubic_x_Root
Cubic_Matrix[1,3]=Cubic_x_Tip
Cubic_Matrix[2,3]=Cubic_x_Root
Cubic_Matrix[3,3]=1
Cubic_Matrix[4,3]=1
Cubic_Matrix[1,4]=1
Cubic_Matrix[2,4]=1
Cubic_Matrix[3,4]=0
Cubic_Matrix[4,4]=0
Cubic_Solution_Matrix=zeros(4,1)
Cubic_Solution_Matrix[1,1]=Cubic_y_Tip
Cubic_Solution_Matrix[2,1]=Cubic_y_Root
Cubic_Solution_Matrix[3,1]=Smoothing_Addition_Tip_Gradient
Cubic_Solution_Matrix[4,1]=Smoothing_Addition_Root_Gradient
Cubic_Constant_ABCD=Cubic_Matrix\Cubic_Solution_Matrix

PlanformCubic=zeros(Nu_v_PLANFORM)
for l=1:m_Symm_IND
if Planform_Builder==1&&Nu_v_PLANFORM_SMOOTHING[l]>0
PlanformCubic[l]=Cubic_Constant_ABCD[1,1]*(Nu_v_PLANFORM[l]^3)+Cubic_Constant_ABCD[2,1]*(Nu_v_PLANFORM[l]^2)+Cubic_Constant_ABCD[3,1]*Nu_v_PLANFORM[l]+Cubic_Constant_ABCD[4,1]
c_v_PLANFORM_ADDITIONAL[l]=PlanformCubic[l]-c_v_PLANFORM[l]
end
end

# PLANFORM RULE NUMBER 3
# The planform must have no sharp angles or discontinuities due to the addition of the minimum structural chord
Smoothing_Addition_Root_Index=Smoothing_Addition_Root_Index2
Smoothing_Addition_Tip_Index=Smoothing_Addition_Tip_Index2
Nu_v_PLANFORM_SMOOTHING=zeros(Nu_v_PLANFORM)+Nu_v_PLANFORM
Nu_v_PLANFORM_SMOOTHING[findin(Nu_v_PLANFORM_SMOOTHING,Nu_v_PLANFORM_SMOOTHING[Nu_v_PLANFORM_SMOOTHING.<=Smoothing_Addition_Root_Index])]=0
Nu_v_PLANFORM_SMOOTHING[findin(Nu_v_PLANFORM_SMOOTHING,Nu_v_PLANFORM_SMOOTHING[Nu_v_PLANFORM_SMOOTHING.>=Smoothing_Addition_Tip_Index])]=0
Smoothing_Addition_Root_IndexLocator=findfirst(Nu_v_PLANFORM_SMOOTHING)
Smoothing_Addition_Tip_IndexLocator=findlast(Nu_v_PLANFORM_SMOOTHING)
Smoothing_Addition_Tip_Gradient=((c_v_PLANFORM[Smoothing_Addition_Root_IndexLocator]+c_v_PLANFORM_ADDITIONAL[Smoothing_Addition_Root_IndexLocator])-(c_v_PLANFORM[Smoothing_Addition_Root_IndexLocator-1]+c_v_PLANFORM_ADDITIONAL[Smoothing_Addition_Root_IndexLocator-1]))/(Nu_v_PLANFORM[Smoothing_Addition_Root_IndexLocator]-Nu_v_PLANFORM[Smoothing_Addition_Root_IndexLocator-1])
Smoothing_Addition_Root_Gradient=((c_v_PLANFORM[Smoothing_Addition_Tip_IndexLocator+1]+c_v_PLANFORM_ADDITIONAL[Smoothing_Addition_Tip_IndexLocator+1])-(c_v_PLANFORM[Smoothing_Addition_Tip_IndexLocator]+c_v_PLANFORM_ADDITIONAL[Smoothing_Addition_Tip_IndexLocator]))/(Nu_v_PLANFORM[Smoothing_Addition_Tip_IndexLocator+1]-Nu_v_PLANFORM[Smoothing_Addition_Tip_IndexLocator])
Cubic_x_Tip=Nu_v_PLANFORM[Smoothing_Addition_Root_IndexLocator]
Cubic_x_Root=Nu_v_PLANFORM[Smoothing_Addition_Tip_IndexLocator]
Cubic_y_Tip=c_v_PLANFORM[Smoothing_Addition_Root_IndexLocator]+c_v_PLANFORM_ADDITIONAL[Smoothing_Addition_Root_IndexLocator]
Cubic_y_Root=c_v_PLANFORM[Smoothing_Addition_Tip_IndexLocator]+c_v_PLANFORM_ADDITIONAL[Smoothing_Addition_Tip_IndexLocator]

Cubic_Matrix=zeros(4,4)
Cubic_Matrix[1,1]=Cubic_x_Tip^3
Cubic_Matrix[2,1]=Cubic_x_Root^3
Cubic_Matrix[3,1]=3*Cubic_x_Tip^2
Cubic_Matrix[4,1]=3*Cubic_x_Root^2
Cubic_Matrix[1,2]=Cubic_x_Tip^2
Cubic_Matrix[2,2]=Cubic_x_Root^2
Cubic_Matrix[3,2]=2*Cubic_x_Tip
Cubic_Matrix[4,2]=2*Cubic_x_Root
Cubic_Matrix[1,3]=Cubic_x_Tip
Cubic_Matrix[2,3]=Cubic_x_Root
Cubic_Matrix[3,3]=1
Cubic_Matrix[4,3]=1
Cubic_Matrix[1,4]=1
Cubic_Matrix[2,4]=1
Cubic_Matrix[3,4]=0
Cubic_Matrix[4,4]=0
Cubic_Solution_Matrix=zeros(4,1)
Cubic_Solution_Matrix[1,1]=Cubic_y_Tip
Cubic_Solution_Matrix[2,1]=Cubic_y_Root
Cubic_Solution_Matrix[3,1]=Smoothing_Addition_Tip_Gradient
Cubic_Solution_Matrix[4,1]=Smoothing_Addition_Root_Gradient
Cubic_Constant_ABCD=Cubic_Matrix\Cubic_Solution_Matrix

PlanformCubic=zeros(Nu_v_PLANFORM)
for l=1:m_Symm_IND
if Planform_Builder==1&&Nu_v_PLANFORM_SMOOTHING[l]>0
PlanformCubic[l]=Cubic_Constant_ABCD[1,1]*(Nu_v_PLANFORM[l]^3)+Cubic_Constant_ABCD[2,1]*(Nu_v_PLANFORM[l]^2)+Cubic_Constant_ABCD[3,1]*Nu_v_PLANFORM[l]+Cubic_Constant_ABCD[4,1]
c_v_PLANFORM_ADDITIONAL[l]=PlanformCubic[l]-c_v_PLANFORM[l]
end
end

# PLANFORM RULE NUMBER 4
# The planform must include the fuselage for holistic aerodynamic analysis
c_v_FUSELAGE=Fuselage_Length_Sizing_Factor*(cos(Nu_v_PLANFORM).^Fuselage_Width_Sizing_Factor)

c_v_PLANFORM_ORIGINAL=zeros(c_v_PLANFORM)+c_v_PLANFORM
c_v_PLANFORM=c_v_PLANFORM+c_v_PLANFORM_ADDITIONAL+c_v_FUSELAGE

# CALLING THE PART 1, PART 2, PART 4 - WASHOUT OPTIMIZER
Nu_Line_GRID,Nu_Line_ACTUAL,Nu_Line_WASHOUT,G_Phi_GRID,G_Phi_ACTUAL,Accuracy_Metric,Nu_v_WASHOUT,c_v_WASHOUT,Gamma_BSLD,G_Phi_WASHOUT,Nu_v_WASHOUT_SCALED,Quarter_Chord_Sweep_Line,LE_Line,TE_Line,Aspect_Ratio,m_Symm_IND,Alpha_WASHOUT,Geometric_Twist_Distribution_Deg_WASHOUT=Washout_Optimizer(FreeStream_Velocity,Sigma_SpanRatio,m_IND,b,c_v_PLANFORM,Quarter_Chord_Sweep_Angle_Deg,Mach_Number,CL_RealWorld,Design_AoA_Deg,AlphaL0_Root_Deg,AlphaL0_Tip_Deg,Geometric_Twist_Distribution_Deg,Plotter_Step_Size,m_Limit_Computer,m_Increment,Accuracy_Increment,Washout_Adjust_Factor_Deg,J2,Minimum_Planform_Chord,Maximum_Planform_Chord,Mass,Gravitational_Constant,Rho,Viscosity)

# END PART 6 - FINAL OPTIMIZATION OF MODIFIED PLANFORM
# END PART 6 - FINAL OPTIMIZATION OF MODIFIED PLANFORM
# END PART 6 - FINAL OPTIMIZATION OF MODIFIED PLANFORM
#
#
#
#
#
#
#
#
#
#
# PART 7 - GRID INDEPENDENCE CONFIRMATION
# PART 7 - GRID INDEPENDENCE CONFIRMATION
# PART 7 - GRID INDEPENDENCE CONFIRMATION

Accuracy_Metric_Singular=maximum(Accuracy_Metric)
println(" ")
println("CHECKING GRID INDEPENDENCE AFTER OPTIMIZATION")
println("Previous Percantage Error= ",Accuracy_Metric_Singular*100," %")
println("m = ",m_IND," for Grid Independence")
println(" ")

#=
# Calling grid independence function
Nu_v_IND,c_v_IND,m_IND_v2,Gna_IND,Nu_Line_GRID,G_Phi_GRID,Nu_Line_ACTUAL,G_Phi_ACTUAL,Alpha_ACTUAL,Geometric_Twist_Distribution_Deg_ACTUAL,Accuracy_Metric=Grid_Independence_Iterator(FreeStream_Velocity,Sigma_SpanRatio,m_IND,b,c_v_WASHOUT,Quarter_Chord_Sweep_Angle_Deg,Mach_Number,CL_RealWorld,Design_AoA_Deg,AlphaL0_Root_Deg,AlphaL0_Tip_Deg,Geometric_Twist_Distribution_Deg_WASHOUT,Plotter_Step_Size,m_Limit_Computer,m_Increment,Accuracy_Increment)
# DOESNT WORK BECAUSE IT CANNOT RUN THE OPTIMIZED PLANFORM AND WASHOUT THROUGH THE m_Limit_Computer VALUE AS THE CORRESPONDING VALUES OF PLANFORM AND WASHOUT DO NOT EXIST YET
# DOESNT WORK BECAUSE IT CANNOT RUN THE OPTIMIZED PLANFORM AND WASHOUT THROUGH THE m_Limit_Computer VALUE AS THE CORRESPONDING VALUES OF PLANFORM AND WASHOUT DO NOT EXIST YET
# DOESNT WORK BECAUSE IT CANNOT RUN THE OPTIMIZED PLANFORM AND WASHOUT THROUGH THE m_Limit_Computer VALUE AS THE CORRESPONDING VALUES OF PLANFORM AND WASHOUT DO NOT EXIST YET

# The way grid independence can be checked is to make all of the above into a function, run the function at m=m_IND as originally assessed, then run the function again at m=m_IND+20 (10 additional points per wing). If the change in results is within tolerance, then grid independence can be said to be achieved, if not, then m=m_IND+20 and rinse and repeat

Accuracy_Metric_Singular=maximum(Accuracy_Metric)
if Accuracy_Metric_Singular<=Accuracy_Increment
println(" ")
println("PASS - GRID INDEPENDENCE REMAINS VALID AFTER OPTIMIZATION")
println("PASS - GRID INDEPENDENCE REMAINS VALID AFTER OPTIMIZATION")
println("PASS - GRID INDEPENDENCE REMAINS VALID AFTER OPTIMIZATION")
println(" ")
println("New Percantage Error= ",Accuracy_Metric_Singular*100," %")
println("m = ",m_IND_v2," for Grid Independence")
println(" ")
else
println(" ")
println("FAIL - GRID INDEPENDENCE IS NO LONGER VALID AFTER OPTIMIZATION")
println("FAIL - GRID INDEPENDENCE IS NO LONGER VALID AFTER OPTIMIZATION")
println("FAIL - GRID INDEPENDENCE IS NO LONGER VALID AFTER OPTIMIZATION")
println(" ")
println("New Percantage Error= ",Accuracy_Metric_Singular*100," %")
println("m = ",m_IND_v2," for Grid Independence")
println(" ")
end
=#

# Find the MAC and AC(y) of the final wing
Total_Area=0
for p=1:m_Symm_IND-1
y_Distance=abs(Nu_v_WASHOUT[p+1]-Nu_v_WASHOUT[p])
x_Distance=(c_v_WASHOUT[p+1]+c_v_WASHOUT[p])/2
Area=y_Distance*x_Distance
Total_Area=Total_Area+Area
end
Total_Area=Total_Area*2
c_y_Squared_Integral=0
for p=1:m_Symm_IND-1
y_Distance=abs(Nu_v_WASHOUT[p+1]-Nu_v_WASHOUT[p])
x_Distance=(((c_v_WASHOUT[p+1]^2)+(c_v_WASHOUT[p]^2))/2)
Area=y_Distance*x_Distance
c_y_Squared_Integral=c_y_Squared_Integral+Area
end
MAC=(2/Total_Area)*c_y_Squared_Integral
# MAC CHECKER
#=c_y_Squared_Integral=0
for p=1:m_Symm_IND-1
y_Distance=abs(Nu_v_WASHOUT[p+1]-Nu_v_WASHOUT[p])
x_Distance=(((c_v_WASHOUT[p+1]^2)+(c_v_WASHOUT[p]^2))/2)
Area=y_Distance*x_Distance
c_y_Squared_Integral=c_y_Squared_Integral+Area
end
c_y_Integral=0
for p=1:m_Symm_IND-1
y_Distance=abs(Nu_v_WASHOUT[p+1]-Nu_v_WASHOUT[p])
x_Distance=((c_v_WASHOUT[p+1]+c_v_WASHOUT[p])/2)
Area=y_Distance*x_Distance
c_y_Integral=c_y_Integral+Area
end
MAC_Checker=c_y_Squared_Integral/c_y_Integral=#

# Calculating Cl_Alpha based on Gamma_BSLD. This should be changed in future to be based on Gna_WASHOUT to ensure that if Gna_WASHOUT does not align with Gamma_BSLD (for example if the planform or washout limit was reached), then the actual Gamma Distribution is used rather than the theoretical Gamma Distribution
Cl_Alpha_y=Gamma_BSLD./c_v_PLANFORM_ORIGINAL

Cm_Alpha_y=(Cm_Alpha_Tip-Cm_Alpha_Root).*Nu_v_WASHOUT+(Cm_Alpha_Tip-(Cm_Alpha_Tip-Cm_Alpha_Root).*Nu_v_WASHOUT[1])

Cm_Cl=Cm_Alpha_y./Cl_Alpha_y
dCm_dCl=zeros(m_Symm_IND-1)
x_AC_Nu_v_Plotter=zeros(m_Symm_IND-1)
for i=1:length(dCm_dCl)
dCm_dCl[i]=(Cm_Cl[i+1]-Cm_Cl[i])/(Nu_v_WASHOUT[i+1]-Nu_v_WASHOUT[i])
x_AC_Nu_v_Plotter[i]=(Nu_v_WASHOUT[i]+Nu_v_WASHOUT[i+1])/2
end

x_AC=-dCm_dCl*MAC
x_AC[findin(x_AC,x_AC[x_AC.>=0])]=0
x_AC[findin(x_AC,x_AC[eachindex(x_AC).<=findin(x_AC,maximum(x_AC))])]=maximum(x_AC)
x_AC=x_AC-maximum(x_AC)
x_AC=cat(1,x_AC[1],x_AC[:])
x_AC=-x_AC*(((Quarter_Chord_Sweep_Line[length(Quarter_Chord_Sweep_Line)-1]-Quarter_Chord_Sweep_Line[length(Quarter_Chord_Sweep_Line)])/(Nu_v_WASHOUT[length(x_AC)-1]-Nu_v_WASHOUT[length(x_AC)]))/((x_AC[length(x_AC)-1]-x_AC[length(x_AC)])/(Nu_v_WASHOUT[length(x_AC)-1]-Nu_v_WASHOUT[length(x_AC)])))
Quarter_Chord_Sweep_Line=Quarter_Chord_Sweep_Line+x_AC-c_v_FUSELAGE*Fuselage_x_Offset_Index
LE_Line=Quarter_Chord_Sweep_Line+0.25*c_v_WASHOUT
TE_Line=Quarter_Chord_Sweep_Line-0.75*c_v_WASHOUT

using Plots
plotly()
plot()
plot!(Nu_v_WASHOUT,c_v_WASHOUT)
plot!(Nu_v_WASHOUT,Gamma_BSLD)
plot!(Nu_Line_WASHOUT,G_Phi_WASHOUT)
gui()

using Plots
plotly()
plot()
xlims!(-1.1*b/2,1.1*b/2)
ylims!(-0.75*1.1*b/2,0.75*1.1*b/2)
plot!(Nu_v_WASHOUT_SCALED,Quarter_Chord_Sweep_Line)
plot!(Nu_v_WASHOUT_SCALED,LE_Line)
plot!(Nu_v_WASHOUT_SCALED,TE_Line)
plot!(-Nu_v_WASHOUT_SCALED,Quarter_Chord_Sweep_Line)
plot!(-Nu_v_WASHOUT_SCALED,LE_Line)
plot!(-Nu_v_WASHOUT_SCALED,TE_Line)
plot!(Nu_v_WASHOUT_SCALED,c_v_WASHOUT+1)
plot!(Nu_v_WASHOUT_SCALED,0.15.*ones(m_Symm_IND,1)+1) # 150mm mark for acceptable Reynolds Number
plot!(-Nu_v_WASHOUT_SCALED,c_v_WASHOUT+1)
plot!(-Nu_v_WASHOUT_SCALED,0.15.*ones(m_Symm_IND,1)+1) # 150mm mark for acceptable Reynolds Number
gui()

using Plots
plotly()
plot()
plot!(Nu_v_WASHOUT,Alpha_WASHOUT*180/pi)
plot!(Nu_v_WASHOUT,Geometric_Twist_Distribution_Deg_WASHOUT)
gui()

using Plots
plotly()
plot()
plot!(Nu_v_PLANFORM,c_v_PLANFORM)
plot!(Nu_v_PLANFORM,c_v_PLANFORM_ORIGINAL)
plot!(Nu_v_PLANFORM,c_v_PLANFORM_ADDITIONAL)
plot!(Nu_v_PLANFORM,c_v_FUSELAGE)
gui()

using Plots
plotly()
plot()
xlims!(0,1)
ylims!(-1,5)
plot!(Nu_v_WASHOUT,Cm_Alpha_y)
plot!(Nu_v_WASHOUT,Cl_Alpha_y)
plot!(Nu_v_WASHOUT,Cm_Cl)
gui()

using Plots
plotly()
plot()
xlims!(0,1)
ylims!(-0.4,0.4)
plot!(Nu_v_WASHOUT,x_AC)
gui()






# ESSENTIAL MODIFICATIONS
# 1. Verify against NACA Report 1208 and the Prandtl-D specifications, then validate specific wing against OpenFOAM simulation
# 2. Calculate the dihedral quarter chord line
# 3. Add in a binary variable for each cubic smoothing line to easily turn them on or off

# DESIRABLE MODIFICATIONS
# 1. Make plots to scale and label them accordingly
# 2. Implement a validation loop of the grid independence study by rerunning the grid independence function with the final wing. If still valid then the solution is final, if not valid, then re-run grid independence study and planform and washout functions on the new m-value attained
# 3. Make a function to rake the tip back
# 4. Double check the weight vs lift calculation by integrating the Gamma curve and convert to lift. This should be equal to the input weight
# 5. Incorporate a changing 3D lift curve slope value across the span in the same way AlphaL0 changes with span location
# 6. Make an automatic output text file

# POSSIBLE BUGS
# 1. Still potential for the while loop to end when residuals are less than 1.0 (occurs in both Planform Optimizer and Washout Optimizer) (seems to be most prevalent when the m-value is very small (less than ~19))
# 2. Discrepency between the calculated Cl_Alpha_y value due to the use of Gamma_BSLD rather than Gna_WASHOUT. This could induce an error if the planform isnt optimized (highly likely if the planform and washout limits are hit, or on extreme sweep angles or aspect ratios)
# 3. Possible error in the way kv is calculated (Beta is eliminated)







#
