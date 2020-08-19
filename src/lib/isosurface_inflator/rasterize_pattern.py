import subprocess, argparse, tempfile, os

parser = argparse.ArgumentParser(description='Rasterize a pattern to a voxel grid')
parser.add_argument('symmetry',    type=str, help='symmetry type (orthotropic, 2d_orthotropic, triply_periodic, 2d_doubly_periodic, etc.)')
parser.add_argument('topology',    type=str, help='path to pattern topology line mesh')
parser.add_argument('densityPath', type=str, help='output path for the density field')

parser.add_argument('--resolution',       type=int,   default=256,    help='resolution of rasterization grid')
parser.add_argument('--voidDensity',      type=float, default=0.0001, help='material density used to void cells')
parser.add_argument('--parameters',       type=str,   default=None,   help='pattern parameters (whitespace separated string)')
parser.add_argument('--triangulatedMesh', type=str,   default=None,   help='output path for the triangulated solid voxels')

args = parser.parse_args();

dim = 2 if (args.symmetry[:2].lower() == "2d") else 3
resString = 'x'.join([str(args.resolution)] * dim)

tmpPath = next(tempfile._get_candidate_names()) + '.msh'

cmd = ['isosurface_cli', args.symmetry, args.topology, '--rasterize', tmpPath, '--rasterResolution', resString]
if (args.parameters): cmd += ['--params',  args.parameters]

subprocess.call(['isosurface_cli', args.symmetry, args.topology, '--rasterize', tmpPath, '--rasterResolution', resString])
subprocess.call(['msh_processor', tmpPath, '--push', str(args.voidDensity), '--push', str(1.0 - args.voidDensity), '-e', 'indicator', '--mul', '--add', '--rename', 'density', '-o', args.densityPath])

if (args.triangulatedMesh):
    subprocess.call(['msh_processor', tmpPath, '-e', 'indicator', '--filterElements', '-o', args.triangulatedMesh])
    subprocess.call(['mesh_convert', args.triangulatedMesh, '-q0', args.triangulatedMesh])

os.remove(tmpPath)
