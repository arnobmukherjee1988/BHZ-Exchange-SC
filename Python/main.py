from utils import *
from config import *
import GenerateSpin as spin
import GenerateHamiltonian as ham
import Observables as obs
import Plotting as plot

def main():
  prog_start_time = time.time()
  
  # Check if the file exists
  if os.path.exists("runtimes.log"):
    os.remove("runtimes.log")
  if os.path.exists("output.log"):
    os.remove("output.log")

  # Create the directory
  os.makedirs('Data', exist_ok=True)
  os.makedirs('Plot', exist_ok=True)

  # Print simulation parameters
  print_parameters(print_info = True)
  
  # Calling, saving, and plotting of the spin configuration. Lastly, calculating skyrmion number
  Sx, Sy, Sz, beta_opt = spin.get_spin()
  datafile = 'spin'
  save_matrices_to_file ('./Data/'+datafile+'.txt', Sx, Sy, Sz)
  plot.SpinPlot(Mat1=Sx, Mat2=Sy, Mat3=Sz, arrow_length=1.5, SaveAs='./Plot/'+datafile+'.pdf')
  Skyrmion_Number = spin.CalculateSkyrmionNumber(Sx, Sy, Sz)
  outputfile = open('output.log', 'a+')
  outputfile.write(f"Skyrmion_Number = {Skyrmion_Number}\n")

  # Creating the Hamiltonian matrix
  start_time = time.time()
  Ham_Matrix = ham.BHZ_hamiltonian(Sx, Sy, Sz)
  # plot.visualize_matrices ([Ham_Matrix], ['Ham_Matrix'], aspect='equal')
  end_time = time.time()
  log_runtime("Creating the Hamiltonian matrix", end_time - start_time)

  # Release spin configuration variables
  del Sx, Sy, Sz
  
  # Diagonalization of the Ham_Matrix to get the eigenvalue array W(n) and eigenfunction matrix Z(n,n)
  start_time = time.time()
  Eigenval, Eigenvec = np.linalg.eigh(Ham_Matrix)
  end_time = time.time()
  log_runtime("Diagonalization of the Hamiltonian matrix", end_time - start_time)
  
  # Release Hamiltonian matrix variable
  del Ham_Matrix
  
  # Save and plot the eigenvalues W in datafile
  datafile = 'eigenvals'
  save_matrices_to_file ('./Data/'+datafile+'.txt', Eigenval)
  # X, Y = np.arange(len(W))[::20], W[::20]
  X, Y = np.arange(len(Eigenval))[:], Eigenval[:]
  plot.plot_data([[X,Y]], xlim=None, ylim=None, linetype='', pointtype='.', 
                  pointsize=20, color='red', xlabel='Eigen Index', ylabel='Eigen Values',
                  title='Eigen Spectrum', SaveAs='./Plot/'+datafile+'.pdf')
  
  # Calculating the Quadrupole Moment
  start_time = time.time()
  QM = obs.QuadrupoleMoment (Eigenvec)
  outputfile.write(f"Quadrupole_Moment = {QM}\n")
  end_time = time.time()
  log_runtime("Calculating the QuadrupoleMoment", end_time - start_time)
  
  # Calculate total DoS
  start_time = time.time()
  Energy, DoS = obs.dos_cal (Eigenval)
  datafile = 'TDoS'
  save_matrices_to_file ('./Data/'+datafile+'.txt', Energy, DoS)
  plot.plot_data([[Energy,DoS]], xlim=None, ylim=None, linetype='-', pointtype='', 
                  pointsize=20, color='red', xlabel='Energy', ylabel='DOS',
                  title='Total DOS', SaveAs='./Plot/'+datafile+'.pdf')
  end_time = time.time()
  log_runtime("Calculating Total DoS", end_time - start_time)
  
  # Calculate LDOS for all sites at the energy of interest and Print
  start_time = time.time()
  LDoS = obs.LDOS_ForAllSites_AtFixedEnergy (energy_of_interest, Eigenval, Eigenvec)
  datafile = 'LDoS'
  save_matrices_to_file ('./Data/'+datafile+'.txt', LDoS)
  plot.visualize_matrices ([LDoS], ['Local DoS'], aspect='equal', interpolation='bicubic', SaveAs = './Plot/'+datafile+'.pdf')
  end_time = time.time()
  log_runtime("Calculating local DoS", end_time - start_time)
  
  del LDoS, Eigenvec, Eigenval

  prog_end_time = time.time()
  log_runtime("Total runtime", prog_end_time - prog_start_time)
  
  outputfile.close()            
  os.system ("mv *.log Data/")

if __name__ == "__main__":
  #start_time = time.time()
  main()
  #end_time = time.time()
  #log_runtime("Total runtime", end_time - start_time)