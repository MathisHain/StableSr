# StableSr
# Inverse inverse/forward model of stable Sr isotope data of Payton and colleagues
# Author: Mathis Hain, mhain at ucsc.edu, mathis-hain.net
# Reference: Paytan et al. (2021), A 35 Myr Record of Seawater Stable Sr Isotopes Reveals a Fluctuating Global Carbon Cycle, Science, doi: TBD

To run the model with the same parameter settings as in Paytan et al (2021) and then generate the figure you can direcly run the "InverseModel_Vfinal.py" by typing:
>> python InverseModel_Vfinal.py

########### USAGE GUIDE #####################
if __name__ == "__main__": # execute the code on the comand line >> python InverseModel_Vfinal.py
    import InverseModel_Vfinal as ModelLib # imports itself
    ModelObject = ModelLib.StableStrontium() # run the model as in Paytan et al 2021
    ModelObject.PublishedFigure() # generate the pulished figure
