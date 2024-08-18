import os, sys

# Set the directory containing the files
directory = sys.argv[1]

i = 0
j = 1
command = "hadd mc_v14-8gev-8.0GeV-1e-ecal_photonuclear_run1.root"
# Loop through each file in the directory
os.chdir(directory)
for filename in os.listdir(directory):
    if "trig" in filename : continue
    if not "histo_v4"  in filename : continue
    i += 1
    # Check if the current file is a file (not a directory)
    if os.path.isfile(os.path.join(directory, filename)):
        # Do something with the file, for example, print its name
        spaceFileName = " " + str(filename)
        command += spaceFileName
    if i == 500 :
        print(command)
        os.system(command)
        i = 0
        j += 1
        command = "hadd mc_v14-8gev-8.0GeV-1e-ecal_photonuclear_run" + str(j)+".root"
