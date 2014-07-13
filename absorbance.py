import numpy as np

# usage: python absorbance.py

background_file=""

input_file="results/C_M_out_w1728_m50000_t4000000-spce-4ns"

write_fft=open('results/C_M_bulk_4ns_fft','w+')
write_abs=open('results/abs_bulk_4ns_fft','w+')
cat_cm_data=open('results/C_M_bulk_4ns_reflected','w+')

data_in=np.loadtxt(input_file,skiprows=0,usecols=(0,1))

cm=data_in[:,1]
plot_time=data_in[:,0]
if background_file:
    print "subtracting background"
    background_in=np.loadtxt(background_file,skiprows=0,usecols=(0,1))
    background_cm=background_in[:,1]
    for i in np.arange(len(cm)):
        cm[i]=background_cm[i]-cm[i]
else:
    print "no background file found"
rev_cm=cm[::-1]

cat_cm=np.zeros(2*len(cm))
cat_time=np.zeros(2*len(cm))
for i in np.arange(len(cm)):
    cat_cm[i]=rev_cm[i]
    cat_cm[i+len(cm)]=cm[i]
    cat_time[i]=-plot_time[len(plot_time)-i-1]
    cat_time[i+len(cm)]=plot_time[i]

transform=np.fft.fft(cat_cm)

n=cat_cm.size
timestep=1.0/1000.0
freq=np.fft.fftfreq(n,d=timestep)

print "frequencies:"
print freq

length=len(cm)
for i in np.arange(length):
    write_fft.write("%f %f\n" % (freq[i],abs(transform[i])))
    write_abs.write("%f %f\n" % (freq[i],freq[i]*freq[i]*abs(transform[i])))
for i in np.arange(2*length):
    cat_cm_data.write("%f %f\n" % (cat_time[i],cat_cm[i]))
