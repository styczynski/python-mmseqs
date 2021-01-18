import os
import pathlib
from typing import Optional, Sequence, Tuple

from .unafold_python_native import HybridMinArgs, _hybrid_min

DEFAULT_MFOLD = [5, "*", 100]


#
# -V, --version
# -h, --help
# -n, --NA=(RNA | DNA) (defaults to RNA)
# -t, --tmin=<minimum temperature> (defaults to 37)
# -i, --tinc=<temperature increment> (defaults to 1)
# -T, --tmax=<maximum temperature> (defaults to 37)
# -s, --suffix=<free energy suffix>
# -o, --output=<prefix>
# -N, --sodium=<[Na+] in M> (defaults to 1)
# -M, --magnesium=<[Mg++] in M> (defaults to 0)
# -p, --polymer
# -r, --prohibit=<i,j,k>
# -E, --energyOnly
# -I, --noisolate
# -z, --zip
# -F, --mfold[=<P,W,MAX>] (defaults to 5,*,100; W determined by sequence length)
# -q, --quiet
# -c, --constraints[=<name of constraint file>] (defaults to prefix.aux)
# -b, --basepairs=<name of basepairs file>
#
def hybrid_min_files(
    files: Sequence[str],
    na: Optional[str] = None,
    tmin: Optional[float] = None,
    tinc: Optional[float] = None,
    tmax: Optional[float] = None,
    suffix: Optional[str] = None,
    output: Optional[str] = None,
    sodium: Optional[float] = None,
    magnesium: Optional[float] = None,
    polymer: bool = False,
    prohibit: Optional[str] = None,
    mfold: Optional[Sequence[Tuple[str, str, str]]] = None,
    constraints: Optional[str] = None,
    basepairs: Optional[str] = None,
    quiet: bool = False,
) -> Tuple[float, float, float]:
    args = HybridMinArgs()
    #
    # if prohibit is None:
    #     prohibit = []
    #
    if na is not None:
        args.na = na
    if tmin is not None:
        args.tmin = tmin
    if tinc is not None:
        args.tinc = tinc
    if tmax is not None:
        args.tmax = tmax
    if suffix is not None:
        args.suffix = suffix
    if output is not None:
        args.output = output
    if sodium is not None:
        args.sodium = sodium
    if magnesium is not None:
        args.magnesium = magnesium
    if polymer:
        args.polymer = True
    if mfold is not None:
        if len(mfold) > 0:
            pstr = ",".join(map(lambda x: str(x), mfold))
            args.mfold.append(pstr)
    if constraints is not None:
        args.constraints = constraints
    if quiet:
        args.quiet = quiet

    args.input1 = files[0]
    args.input2 = files[1]

    data_path = os.path.abspath(
        os.path.join(pathlib.Path(__file__).parent.absolute(), "..", "data")
    )
    if not os.path.exists(data_path):
        data_path = os.path.abspath(
            os.path.join(
                pathlib.Path(__file__).parent.absolute(), "..", "..", "data"
            )
        )
        if not os.path.exists(data_path):
            raise Exception(
                "Failed to determine UNAFOLD data directory."
                f"Checked from: {pathlib.Path(__file__).parent.absolute()}"
            )
    os.environ["UNAFOLDDAT"] = data_path

    ret = _hybrid_min(args)
    if ret[0] != 0.0:
        raise Exception("Failed to calculate hybrid_min.")
    return ret[1:]


def hybrid_min(
    sequences: Sequence[str],
    na: Optional[str] = None,
    tmin: Optional[float] = None,
    tinc: Optional[float] = None,
    tmax: Optional[float] = None,
    suffix: Optional[str] = None,
    output: Optional[str] = None,
    sodium: Optional[float] = None,
    magnesium: Optional[float] = None,
    polymer: bool = False,
    prohibit: Optional[str] = None,
    mfold: Optional[Sequence[Tuple[str, str, str]]] = None,
    constraints: Optional[str] = None,
    basepairs: Optional[str] = None,
) -> Tuple[float, float, float]:
    return hybrid_min_files(
        sequences,
        na=na,
        tmin=tmin,
        tinc=tinc,
        tmax=tmax,
        suffix=suffix,
        output=output,
        sodium=sodium,
        magnesium=magnesium,
        polymer=polymer,
        prohibit=prohibit,
        mfold=mfold,
        constraints=constraints,
        basepairs=basepairs,
        quiet=True,
    )
