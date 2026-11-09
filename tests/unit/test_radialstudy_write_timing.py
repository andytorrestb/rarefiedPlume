from pathlib import Path
import re
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[2]/'tests/unit'))
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from test_radialstudy_physics import model
from plumetools.radialstudy.case import write_case
from plumetools.radialstudy.mesh import initial_mesh

def test_native_write_interval(tmp_path):
    m = model(100)
    case = tmp_path/'case'
    record = write_case(case, m, initial_mesh(m))
    control = (case/'system/controlDict').read_text()
    mode = re.search(r'writeControl\s+(\w+);', control)[1]
    interval = float(re.search(r'writeInterval\s+([^;]+);', control)[1])
    dt = record['numerics']['delta_t_s']
    assert mode == 'runTime', 'hyStrath requires a physical write interval'
    assert round(interval/dt) == record['numerics']['steps_per_window']
    assert interval/dt < 2**31-1, 'native 32-bit label must not overflow'
