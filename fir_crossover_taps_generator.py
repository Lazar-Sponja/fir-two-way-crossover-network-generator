import numpy as np 
import scipy.signal as signal
import argparse, os

# Tap estimator algorithms
def harris_approx(Bt, Aa, fs=None) -> int:
    """Function for approximating FIR order using Harris's rule of thumb

    Args:
        Wp (float): Band transition width. If fs is not None, value in Hz is expected.
        Aa (float): Desired stopband attenuation.
        fs (float, optional): Sampling frequency of digital system. Defaults to None.

    Returns:
        int: Approximate number of taps needed to satisfy attenuation spec
    """
    Bt /= fs if fs is not None else Bt
    return int(np.ceil(Aa/(22*Bt)))

def bellanger_approx(Bt, dp, da, fs=None) -> int:
    """

    Args:
        Wp (float): Band transition width If fs is not None, value in Hz is expected.
        dp (float): Passband ripple.
        da (float): Stopdand ripple.
        fs (float, optional): Sampling frequency of digital system. Defaults to None.

    Returns:
        int: Approximate number of taps needed to satisfy ripple spcs
    """
    Bt /= fs if fs is not None else Bt
    return int(np.ceil(2*np.log10(1/(10*dp*da))/(3*Bt)))

# Arugmet parser instance
parser = argparse.ArgumentParser(description='Generates a text file with taps for 2 way crossover network made of FIR filters')

# parser arguments
exclusive_group1 = parser.add_mutually_exclusive_group(required=True)
exclusive_group2 = parser.add_mutually_exclusive_group(required=True)
exclusive_group1.add_argument(
    '-fc', '--crossover_frequency', 
    help='Crossover frequency. If -fs is not set, normalization by sampling frequency is assumed. This argument is mandatory when -fp and -fa are not set.', 
    type=float,
)
exclusive_group2.add_argument(
    '-Bt', '--transition_width', 
    help='Width of frequency band transition. If -fs is not set, normalization by sampling frequency is assumed. This argument is mandatory when -fp and -fa are not set.', 
    type=float,
)
exclusive_group1.add_argument(
    '-fp', '--passband_frequency', 
    help='Passband frequency of filter. If -fs is not set, normalization by sampling frequency is assumed. This argument is mandatory when -fc and -Bt are not set.', 
    type=float,
)
exclusive_group2.add_argument(
    '-fa', '--stopband_frequency', 
    help='Width of frequency band transition. If fs is not given, normalization by sampling frequency is assumed', 
    type=float,
)
parser.add_argument(
    '-Ap', '--passband-attenuation', 
    help='Attenuation permissible in the passband. Value is in dB', 
    type=float
)
parser.add_argument(
    '-Aa', '--stopband-attenuation', 
    help='Attenuation of the stopband. Value is in dB', 
    type=float
)
parser.add_argument(
    '-N', '--numtaps', 
    help='Number of taps to create filter with. Will override numtap finding algorithm choice', 
    type=float
)
parser.add_argument(
    '--numtaps-finder', 
    help='Algorithm to be used to determine numtaps. Availale algoritms are Harris\'s rule of thumb and Bellanger\'s rule of thumb from "Digital Processing of Signals – Theory and Practice". Defaults to harris if Ap is not give, else bellanger.', 
    type=str, 
    choices=['harris', 'bellanger']
)
parser.add_argument(
    '--fir-algorithm', 
    help='Algorithm used to generate taps. Options are firls and remez. Defaults to firls', 
    type=str, 
    default='firls',
    choices=['firls', 'remez']
)
parser.add_argument(
    '-w', '--weigh-taps', 
    help='Apply weighing to numtaps calculations', 
    action='store_true'
)
parser.add_argument(
    '-o', '--output_path', 
    help='Set output path, optional default to cwd', 
    type=str
)
parser.add_argument(
    '-fs', '--sampling-frequency', 
    help='Sampling frequency. If not specified, inputs are assumed to be normalized by sampling frequency. If specified, input frequencies need to be in Hz.', 
    default=None, 
    type=float
)
args = parser.parse_args()

# Parse output path
if args.output_path is not None: 
    if os.path.isdir(args.output_path):
        lowpass = os.path.join(args.output_path, 'lowpass_crossover.txt')
        highpass = os.path.join(args.output_path, 'highpass_crossover.txt')
else:
    lowpass = os.path.join(os.getcwd(), 'lowpass_crossover.txt')
    highpass = os.path.join(os.getcwd(), 'highpass_crossover.txt')

#Parse if either N or Aa are given
if args.numtaps is None and args.stopband_attenuation is None:
    raise ValueError('Either --numtaps or --stopband-attenuation must be provided')

#Parsing numtaps finder 
if args.numtaps_finder is None:
    numtaps_finder = 'harris' if args.passband_attenuation is None else 'bellanger'
else:
    numtaps_finder = args.numtaps_finder

# Defining passband and stopband frequencies
if args.crossover_frequency is not None:
    Bt = args.transition_width
    fp = args.crossover_frequency - 0.5*Bt
    fa = args.crossover_frequency + 0.5*Bt
else:
    fp = args.passband_frequency
    fa = args.stopband_frequency
    Bt = fa-fp

#Calculating numtaps
if args.numtaps is None:
    match numtaps_finder:
        case 'harris':
            N  = harris_approx(Bt, args.stopband_attenuation, args.sampling_frequency)
        case 'bellanger':
            if args.passband_attenuation:
                dp = (10**(0.05*args.passband_attenuation)-1)/(10**(0.05*args.passband_attenuation)+1) #passband ripple
                da = 10**(-0.05*args.stopband_attenuation) #stopband ripple
                N = bellanger_approx(Bt, dp, da, fs=args.sampling_frequency)
            else:
                raise ValueError('Passband attenuation must be given when using bellanger')
    print(f'numtaps is estimated as {N}')
else:
    N = int(args.numtaps) 

#Calculate weighing values
if args.weigh_taps:
    if args.stopband_attenuation is None or args.passband_attenuation is None:
        raise ValueError('Both stopband and passband attenuations must be given when weighing errors')
    else:
        dp = (10**(0.05*args.passband_attenuation)-1)/(10**(0.05*args.passband_attenuation)+1) #passband ripple
        da = 10**(-0.05*args.stopband_attenuation) #stopband ripple
        W_lp = [1, dp/da]
        W_hp = [dp/da, 1]
else:
    W_lp = None
    W_hp = None

#Calculating taps
F = [0, fp, fa, 0.5*args.sampling_frequency] if args.sampling_frequency else [0, fp, fa, 0.5]
Hid_lp = [1, 1, 0, 0]   #Ideal filter response lp
Hid_hp = [0, 0, 1, 1]   #Ideal filter response hp
match args.fir_algorithm:
    case 'firls':
        h_lp = signal.firls(N if N % 2 else N-1, F, Hid_lp, W_lp, fs=args.sampling_frequency)
        h_hp= signal.firls(N if N % 2 else N-1, F, Hid_hp, W_hp, fs=args.sampling_frequency)
    case 'remez':
        h_lp = signal.remez(N, F, Hid_lp[::2], W_lp, fs=args.sampling_frequency)
        h_hp= signal.remez(N if N % 2 else N-1, F, Hid_hp[::2], W_hp, fs=args.sampling_frequency)
        
#Save taps to txt file
np.savetxt(lowpass, h_lp)
np.savetxt(highpass,h_hp)


