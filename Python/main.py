from utils import *
from config import *
import GenerateSpinNew as spin
import GenerateSpinConfig as spin_old
import GenerateHamiltonian as ham
import Plotting as plot

def main():
  # Print simulation parameters
  print_parameters(print_info = True)    
  
  # Getting Skyrmion spin configuration and calculating skyrmion number
  Sx, Sy, Sz, Theta, Phi, beta_opt = spin_old.skyrmion_optization()
  Skyrmion_Number = spin.CalculateSkyrmionNumber(Sx, Sy, Sz)
  print("Skyrmion_Number = ", Skyrmion_Number)
  
  # # Getting Skyrmion spin configuration and calculating skyrmion number
  # Sx, Sy, Sz, Theta, Phi, beta_opt = spin.get_spin()
  # Skyrmion_Number = spin.CalculateSkyrmionNumber(Sx, Sy, Sz)
  # print("Skyrmion_Number = ", Skyrmion_Number)
  
  
  # saving the spin configuration in datafile
  datafile = SpinConfig_Type+'_SpinConfig.txt'
  save_matrices_to_file (datafile, Sx, Sy, Sz)

  # Plot the spin configuration
  plot.SpinPlot(Mat1=Sx, Mat2=Sy, Mat3=Sz, arrow_length=1.5, SaveAs=SpinConfig_Type+'.pdf')

  # Creating the Hamiltonian matrix
  # Ham_Matrix = ham.hamiltonian(Sx, Sy, Sz)
  # Ham_Matrix = ham.nambu_hamiltonian(Sx, Sy, Sz)
  # plot.visualize_matrices ([Ham_Matrix], ['Ham_Matrix'], aspect='equal')
  # PrintMatrix (Ham_Matrix)

  # diagonalization of the Ham_Matrix to get the eigenvalue array W(n) and eigenfunction matrix Z(n,n)
  # W, Z = np.linalg.eigh(Ham_Matrix)
  
  # QM = ham.QuadrupoleMoment (Z)
  # print(QM)

  # Bott = ham.Bott_Index(Z)
  # print(Bott)
  
  # save the eigenvalues W in datafile
  # datafile = SpinConfig_Type + '_eigenvals.txt'
  # save_matrices_to_file (datafile, W)

  # plot the eigen spectrum (every 20 points)
  # X, Y = np.arange(len(W))[::20], W[::20]
  # X, Y = np.arange(len(W))[:], W[:]
  # plot.plot_data([[X,Y]], xlim=None, ylim=None, linetype='', pointtype='.', 
  #                 pointsize=20, color='red', xlabel='Eigen Index', ylabel='Eigen Values',
  #                 title='Eigen Spectrum', SaveAs=SpinConfig_Type+'_eigenvals.pdf')
  

if __name__ == "__main__":
  main()