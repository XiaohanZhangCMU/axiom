from solve_pde import solve_pde
from utils import *
import sys
from run_dynamics import run_fine_lorenz, run_coarse_lorenz
from pde_model import Lorenz_Model
import argparse


def main(args):

    parser = argparse.ArgumentParser(description='Process task and param file.')
    parser.add_argument('--Train', '-t','-T', help='train or test', type=int, default=1)
    parser.add_argument('--param', '-p','-P', help='parameter file', type=str, default='lorenz_params.txt')
    parser.add_argument('--restore', '-r','-R', help='restore saved model weights', type=bool, default=False)
    args = parser.parse_args()
    configs = load_params(args.param)

    if args.restore:
        configs['restore'] = True
    else:
        configs['restore'] = False

    print(configs)
    if args.Train:
        model = Lorenz_Model(**configs)
        result = solve_pde(model, **configs)
    else:
        fine_traj = run_fine_lorenz(icpt=[-10,5,23], N=1000, **configs)
        coarse_traj = run_coarse_lorenz(icpt = [-10,23], N=1000, **configs)

        np.save('fine_traj.npy',fine_traj)
        np.save('coarse_traj.npy',coarse_traj)

if __name__ == '__main__':
    main(sys.argv)
