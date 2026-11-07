"""Structure preparation: input parsing, 3D embedding, ion/charge variants.

This is the deterministic replacement for the "draw it in GaussView and hit
Clean" step:

* RDKit ETKDGv3 + MMFF94s generates a chemically sensible 3D conformer
  (the same role GaussView's structure builder / Clean button plays).
* Ions / counter-ions are placed at a physically motivated coordination
  distance and relaxed with a force field.
* The final safety net is a cheap Gaussian pre-optimisation
  (config ``chemistry.preopt``, default HF/3-21G opt) before the production
  level of theory, which is what actually prevents most SCF/link blow-ups.

GaussView itself stays a *human review* terminal (``gflow preview``); it is
not required for the automated path.
"""
from __future__ import annotations

import math
import re
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors
from rdkit.Geometry import Point3D

_PT = Chem.GetPeriodicTable()

# ions commonly used for "add an ion" iteration; users may pass any SMILES
COMMON_IONS = {
    "Na+": "[Na+]", "K+": "[K+]", "Li+": "[Li+]", "Cs+": "[Cs+]",
    "Mg2+": "[Mg+2]", "Ca2+": "[Ca+2]", "Zn2+": "[Zn+2]",
    "Cl-": "[Cl-]", "Br-": "[Br-]", "I-": "[I-]", "F-": "[F-]",
    "NH4+": "[NH4+]", "NO3-": "O=[N+]([O-])[O-]",
    "HCO3-": "OC(=O)[O-]", "H3O+": "[OH3+]",
}


# ---------------------------------------------------------------------------
# small helpers
# ---------------------------------------------------------------------------
def element_symbol(z: int) -> str:
    return _PT.GetElementSymbol(int(z))


def xyz_block(coords: list[tuple[str, float, float, float]], comment: str = "") -> str:
    lines = [f"{len(coords)}", comment or "gaussflow structure"]
    lines += [f"{s:2s} {x:14.6f} {y:14.6f} {z:14.6f}" for s, x, y, z in coords]
    return "\n".join(lines) + "\n"


def coords_from_mol(mol: Chem.Mol) -> list[tuple[str, float, float, float]]:
    conf = mol.GetConformer()
    out = []
    for i, atom in enumerate(mol.GetAtoms()):
        p = conf.GetAtomPosition(i)
        out.append((atom.GetSymbol(), float(p.x), float(p.y), float(p.z)))
    return out


def mol_from_coords(coords: list[tuple[str, float, float, float]],
                    charge: int = 0, bonds: bool = False) -> Chem.Mol:
    """Build an RDKit molecule from element/coordinate triples."""
    xyz = f"{len(coords)}\n\n" + "\n".join(
        f"{s} {x:.6f} {y:.6f} {z:.6f}" for s, x, y, z in coords) + "\n"
    mol = Chem.MolFromXYZBlock(xyz)
    if mol is None:
        raise ValueError("无法从坐标构建分子")
    if bonds and mol.GetNumAtoms() > 1:
        try:
            from rdkit.Chem import rdDetermineBonds
            rdDetermineBonds.DetermineBonds(mol, charge=charge)
        except Exception:
            pass
    return mol


def formula_of_coords(coords: list[tuple[str, float, float, float]]) -> str:
    counts: dict[str, int] = {}
    for s, *_ in coords:
        counts[s] = counts.get(s, 0) + 1
    # simple Hill ordering: C, H, then alphabetical
    order = []
    if "C" in counts:
        order.append("C")
    if "H" in counts:
        order.append("H")
    order += sorted(k for k in counts if k not in ("C", "H"))
    parts = []
    for k in order:
        n = counts[k]
        parts.append(k if n == 1 else f"{k}{n}")
    return "".join(parts)


# ---------------------------------------------------------------------------
# 3D embedding
# ---------------------------------------------------------------------------
def embed_3d(mol: Chem.Mol, seed: int = 42, n_conformers: int = 1,
             optimize: bool = True) -> Chem.Mol:
    """Return a molecule with a 3D conformer (ETKDGv3 + force-field cleanup)."""
    mol = Chem.AddHs(Chem.Mol(mol))
    n_conf = max(1, int(n_conformers))
    params = AllChem.ETKDGv3()
    params.randomSeed = int(seed)
    params.useSmallRingTorsions = True
    params.pruneRmsThresh = 0.25 if n_conf > 1 else -1.0
    conf_ids = list(AllChem.EmbedMultipleConfs(mol, numConfs=n_conf, params=params))
    if not conf_ids:
        params.useRandomCoords = True
        params.maxAttempts = 500
        conf_ids = list(AllChem.EmbedMultipleConfs(mol, numConfs=n_conf, params=params))
    if not conf_ids:
        raise ValueError("RDKit 无法为该分子生成 3D 结构（可改用 --file 提供坐标）")

    if optimize:
        ff_name = _optimize_all(mol, conf_ids)
        # pick lowest-energy conformer and make it conformer 0
        best, best_e = conf_ids[0], float("inf")
        for cid in conf_ids:
            try:
                e = _ff_energy(mol, cid, ff_name)
            except Exception:
                e = float("inf")
            if e < best_e:
                best, best_e = cid, e
        if best != conf_ids[0]:
            conf = _conf_copy(mol, best)
            mol.RemoveAllConformers()
            mol.AddConformer(conf, assignId=True)
    return mol


def _ff_energy(mol: Chem.Mol, conf_id: int, ff_name: str) -> float:
    if ff_name == "MMFF":
        props = AllChem.MMFFGetMoleculeProperties(mol, mmffVariant="MMFF94s")
        if props is None:
            raise ValueError("MMFF 参数缺失")
        ff = AllChem.MMFFGetMoleculeForceField(mol, props, confId=conf_id)
    else:
        ff = AllChem.UFFGetMoleculeForceField(mol, confId=conf_id)
    if ff is None:
        raise ValueError("无法建立力场")
    return ff.CalcEnergy()


def _optimize_all(mol: Chem.Mol, conf_ids: list[int]) -> str:
    try:
        props = AllChem.MMFFGetMoleculeProperties(mol, mmffVariant="MMFF94s")
        if props is not None:
            for cid in conf_ids:
                ff = AllChem.MMFFGetMoleculeForceField(mol, props, confId=cid)
                if ff is None:
                    raise ValueError("MMFF force field failed")
                ff.Minimize(maxIts=2000)
            return "MMFF"
    except Exception:
        pass
    try:
        for cid in conf_ids:
            ff = AllChem.UFFGetMoleculeForceField(mol, confId=cid)
            if ff is None:
                raise ValueError("UFF force field failed")
            ff.Minimize(maxIts=2000)
        return "UFF"
    except Exception:
        return "none"


def _conf_copy(mol: Chem.Mol, conf_id: int) -> Chem.Conformer:
    src = mol.GetConformer(conf_id)
    conf = Chem.Conformer(mol.GetNumAtoms())
    for i in range(mol.GetNumAtoms()):
        conf.SetAtomPosition(i, src.GetAtomPosition(i))
    return conf


# ---------------------------------------------------------------------------
# input parsing
# ---------------------------------------------------------------------------
def molecule_from_smiles(smiles: str, seed: int = 42, n_conformers: int = 1) -> Chem.Mol:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"无法解析 SMILES: {smiles}")
    return embed_3d(mol, seed=seed, n_conformers=n_conformers)


def molecule_from_file(path: str | Path, seed: int = 42,
                       charge: int | None = None) -> tuple[Chem.Mol, dict]:
    path = Path(path)
    suffix = path.suffix.lower()
    meta: dict = {"source": str(path)}
    if suffix in (".sdf", ".mol"):
        supplier = Chem.SDMolSupplier(str(path), removeHs=False)
        mol = next((m for m in supplier if m is not None), None)
        if mol is None:
            raise ValueError(f"无法读取结构文件: {path}")
        if mol.GetNumConformers() == 0:
            mol = embed_3d(mol, seed=seed)
        return mol, meta
    if suffix in (".xyz",):
        coords = read_xyz(path)
        mol = mol_from_coords(coords, charge=charge or 0, bonds=True)
        meta["charge"] = charge
        return mol, meta
    if suffix in (".gjf", ".com"):
        parsed = read_gjf(path)
        mol = mol_from_coords(parsed["coords"], charge=parsed["charge"], bonds=True)
        meta.update(charge=parsed["charge"], multiplicity=parsed["multiplicity"],
                    route=parsed["route"], title=parsed["title"])
        return mol, meta
    raise ValueError(f"不支持的结构文件类型: {suffix}")


def read_xyz(path: str | Path) -> list[tuple[str, float, float, float]]:
    lines = Path(path).read_text(encoding="utf-8", errors="replace").splitlines()
    if not lines:
        raise ValueError(f"空文件: {path}")
    n = int(lines[0].split()[0])
    coords = []
    for line in lines[2:2 + n]:
        parts = line.split()
        if len(parts) < 4:
            continue
        coords.append((parts[0].capitalize(), float(parts[1]), float(parts[2]), float(parts[3])))
    if len(coords) != n:
        raise ValueError(f"XYZ 原子数不符: 期望 {n}, 读取 {len(coords)}")
    return coords


def read_gjf(path: str | Path) -> dict:
    """Minimal Gaussian input reader (route, title, charge/mult, coords)."""
    text = Path(path).read_text(encoding="utf-8", errors="replace")
    lines = text.splitlines()
    route = ""
    title = ""
    charge = 0
    mult = 1
    coords: list[tuple[str, float, float, float]] = []
    state = "route"
    coord_lines: list[str] = []
    for line in lines:
        s = line.strip()
        if state == "route":
            if s.startswith("#"):
                route = s.lstrip("#").strip()
                route = re.sub(r"^[pPnNtT]\s+", "", route).strip()
                state = "after_route"
            elif s.startswith("%"):
                continue
            continue
        if state == "after_route":
            if s == "":
                state = "title"
            continue
        if state == "title":
            if s == "":
                state = "charge"
            elif not title:
                title = s
            else:
                title += " " + s
            continue
        if state == "charge":
            m = re.match(r"^(-?\d+)\s+(\d+)$", s)
            if m:
                charge, mult = int(m.group(1)), int(m.group(2))
                state = "coords"
            continue
        if state == "coords":
            if s == "":
                if coord_lines:
                    state = "done"
                continue
            m = re.match(r"^([A-Za-z]{1,2})\s+(-?\d+(?:\.\d*)?(?:[EeDd][-+]?\d+)?)\s+"
                         r"(-?\d+(?:\.\d*)?(?:[EeDd][-+]?\d+)?)\s+"
                         r"(-?\d+(?:\.\d*)?(?:[EeDd][-+]?\d+)?)", s)
            if m:
                def f(x: str) -> float:
                    return float(x.replace("D", "E").replace("d", "e"))
                coord_lines.append(s)
                coords.append((m.group(1).capitalize(), f(m.group(2)), f(m.group(3)), f(m.group(4))))
    if not route:
        raise ValueError(f"未找到 Gaussian 路由行: {path}")
    if not coords:
        raise ValueError(f"未找到坐标: {path}")
    return {"route": route, "title": title, "charge": charge,
            "multiplicity": mult, "coords": coords}


# ---------------------------------------------------------------------------
# charge / multiplicity inference
# ---------------------------------------------------------------------------
def infer_charge_mult(mol: Chem.Mol) -> tuple[int, int]:
    charge = int(Chem.GetFormalCharge(mol))
    electrons = sum(a.GetAtomicNum() for a in mol.GetAtoms()) - charge
    radicals = sum(a.GetNumRadicalElectrons() for a in mol.GetAtoms())
    mult = (radicals + 1) if radicals > 0 else (1 if electrons % 2 == 0 else 2)
    # guard: multiplicity parity must match electron count
    if (electrons - (mult - 1)) % 2 != 0:
        mult = 1 if electrons % 2 == 0 else 2
    return charge, mult


# ---------------------------------------------------------------------------
# ions / variants
# ---------------------------------------------------------------------------
def _anchor_atom(mol: Chem.Mol) -> int:
    """Pick the best coordination site: most polarisable hetero atom."""
    best, best_score = 0, -1.0
    for atom in mol.GetAtoms():
        z = atom.GetAtomicNum()
        score = {8: 3.0, 7: 2.5, 16: 2.0, 15: 1.5, 9: 2.0}.get(z, 0.0)
        score -= 0.4 * atom.GetDegree()
        if score > best_score:
            best, best_score = atom.GetIdx(), score
    return best


def _place(conf: Chem.Conformer, anchor: int, mol: Chem.Mol, charge_ion: int,
           dist: float = 2.35) -> tuple[float, float, float]:
    apos = conf.GetAtomPosition(anchor)
    ax, ay, az = apos.x, apos.y, apos.z
    # outward direction from the heavy-atom neighbours (lone-pair side)
    neigh = [n.GetIdx() for n in mol.GetAtomWithIdx(anchor).GetNeighbors()
             if n.GetAtomicNum() > 1]
    if neigh:
        bx = sum(conf.GetAtomPosition(i).x for i in neigh) / len(neigh)
        by = sum(conf.GetAtomPosition(i).y for i in neigh) / len(neigh)
        bz = sum(conf.GetAtomPosition(i).z for i in neigh) / len(neigh)
        dx, dy, dz = ax - bx, ay - by, az - bz
    else:
        dx, dy, dz = ax, ay, az
    norm = math.sqrt(dx * dx + dy * dy + dz * dz) or 1.0
    dx, dy, dz = dx / norm, dy / norm, dz / norm
    if charge_ion < 0:
        # anion: point away from the lone-pair side, closer to the H-bond donor
        hneigh = [n.GetIdx() for n in mol.GetAtomWithIdx(anchor).GetNeighbors()
                  if n.GetAtomicNum() == 1]
        if hneigh:
            hpos = conf.GetAtomPosition(hneigh[0])
            hx, hy, hz = hpos.x - ax, hpos.y - ay, hpos.z - az
            hn = math.sqrt(hx * hx + hy * hy + hz * hz) or 1.0
            dx, dy, dz = hx / hn, hy / hn, hz / hn
    return ax + dx * dist, ay + dy * dist, az + dz * dist


def add_ion(mol: Chem.Mol, ion: str, seed: int = 42,
            dist: float = 2.35, anchor_index: int | None = None) -> Chem.Mol:
    """Add one ion (SMILES or common name) to a 3D molecule.

    ``anchor_index`` pins the coordination atom (e.g. a chosen O or F for Li⁺);
    otherwise the most polarisable hetero atom is used.
    """
    ion_smiles = COMMON_IONS.get(ion, ion)
    ion_mol = Chem.MolFromSmiles(ion_smiles)
    if ion_mol is None:
        raise ValueError(f"无法解析离子: {ion}")
    ion_charge = int(Chem.GetFormalCharge(ion_mol))
    ion_mol = Chem.AddHs(ion_mol)
    # embed the (polyatomic) ion; single atoms get a trivial conformer
    if ion_mol.GetNumAtoms() > 1:
        params = AllChem.ETKDGv3()
        params.randomSeed = int(seed)
        if AllChem.EmbedMolecule(ion_mol, params) != 0:
            params.useRandomCoords = True
            AllChem.EmbedMolecule(ion_mol, params)
        if ion_mol.GetNumConformers() == 0:
            raise ValueError(f"无法为离子生成 3D 结构: {ion}")
    else:
        cf = Chem.Conformer(1)
        cf.SetAtomPosition(0, Point3D(0.0, 0.0, 0.0))
        ion_mol.RemoveAllConformers()
        ion_mol.AddConformer(cf, assignId=True)

    base_n = mol.GetNumAtoms()
    src = mol.GetConformer()
    anchor = int(anchor_index) if anchor_index is not None else _anchor_atom(mol)
    pos = _place(src, anchor, mol, ion_charge, dist=dist)
    iconf = ion_mol.GetConformer()
    first = iconf.GetAtomPosition(0)
    conf = Chem.Conformer(base_n + ion_mol.GetNumAtoms())
    for i in range(base_n):
        conf.SetAtomPosition(i, src.GetAtomPosition(i))
    for j in range(ion_mol.GetNumAtoms()):
        p = iconf.GetAtomPosition(j)
        conf.SetAtomPosition(base_n + j, Point3D(pos[0] + p.x - first.x,
                                                pos[1] + p.y - first.y,
                                                pos[2] + p.z - first.z))
    combined = Chem.CombineMols(mol, ion_mol)
    combined.RemoveAllConformers()
    combined.AddConformer(conf, assignId=True)
    try:
        Chem.SanitizeMol(combined)  # initialise ring info / property cache
    except Exception:
        pass

    # force-field relaxation (MMFF often lacks alkali/alkaline-earth params)
    optimized = False
    try:
        props = AllChem.MMFFGetMoleculeProperties(combined, mmffVariant="MMFF94s")
        if props is not None:
            ff = AllChem.MMFFGetMoleculeForceField(combined, props, confId=0)
            if ff is not None:
                ff.Minimize(maxIts=1000)
                optimized = True
    except Exception:
        optimized = False
    if not optimized:
        try:
            ff = AllChem.UFFGetMoleculeForceField(combined, confId=0)
            if ff is not None:
                ff.Minimize(maxIts=1000)
        except Exception:
            pass
    return combined


def build_variant(base_mol: Chem.Mol, ions: list[str], seed: int = 42,
                  variant_index: int = 0) -> tuple[Chem.Mol, dict]:
    """Return a (possibly ion-doped) molecule plus variant metadata."""
    if not ions:
        return Chem.Mol(base_mol), {"type": "base", "ions": []}
    ion = ions[variant_index % len(ions)]
    mol = add_ion(base_mol, ion, seed=seed + variant_index)
    return mol, {"type": "ion", "ion": ion, "ions": [ion]}


def list_coordination_sites(mol: Chem.Mol,
                            elements: list[str] | None = None) -> list[dict]:
    """Candidate Li⁺ coordination sites: O (sp²/sp³) and F atoms."""
    elements = [e.upper() for e in (elements or ["O", "F"])]
    sites = []
    for atom in mol.GetAtoms():
        sym = atom.GetSymbol().upper()
        if sym not in elements:
            continue
        # O with at least one heavy neighbour (any O in organics qualifies);
        # F is only interesting when bound to carbon (not ionic F⁻).
        if sym == "F":
            neigh = [n for n in atom.GetNeighbors() if n.GetAtomicNum() > 1]
            if not neigh or neigh[0].GetSymbol().upper() not in ("C", "S", "P", "B"):
                continue
        sites.append({
            "atom_index": atom.GetIdx(),
            "element": sym,
            "degree": atom.GetDegree(),
            "idx_in_element": sum(1 for s in sites if s["element"] == sym),
        })
    return sites
