



Printer_x=zeros(151)
Printer_z=zeros(151)
Printer_xz=zeros(2,301)

for i=1:301

Printer_x_Upper=AIRFOIL_OUTPUT_x_Upper[i,:]
Printer_x_Lower=flipdim(AIRFOIL_OUTPUT_x_Lower[i,:],1)
deleteat!(Printer_x_Lower,151)
Printer_x=cat(1,Printer_x_Lower,Printer_x_Upper)
Printer_x=Printer_x.'

Printer_z_Upper=AIRFOIL_OUTPUT_z_Upper[i,:]
Printer_z_Lower=flipdim(AIRFOIL_OUTPUT_z_Lower[i,:],1)
deleteat!(Printer_z_Lower,151)
Printer_z=cat(1,Printer_z_Lower,Printer_z_Upper)
Printer_z=Printer_z.'

Printer_xz[1,:]=Printer_x
Printer_xz[2,:]=Printer_z

Outfile_xz="RTJ M1 Data xxx zzz Station $i.dat"
Outfile_Opener_xz=open(Outfile_xz,"w")

for j=1:301
Data_Row_x=Printer_xz[1,j]
Data_Row_z=Printer_xz[2,j]
write(Outfile_Opener_xz,"$Data_Row_x\t$Data_Row_z\r","\n")
end

close(Outfile_Opener_xz)
end

Printer_y=zeros(301)
Printer_y=FINAL_y_locations.'.'
Outfile_y="RTJ M1 Data yyy Station All.dat"
Outfile_Opener_y=open(Outfile_y,"w")
print(Outfile_Opener_y,Printer_y)
close(Outfile_Opener_y)
