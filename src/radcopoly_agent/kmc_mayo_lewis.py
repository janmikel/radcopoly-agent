from dataclasses import dataclass
from collections import Counter
import random


@dataclass
class KMCResult:
    chains: list[str]
    mn: float
    mw: float
    dispersity: float
    mean_dp: float
    fraction_m1: float
    fraction_m2: float
    avg_block_m1: float
    avg_block_m2: float


def choose_next_monomer(last, f1, r1, r2):
    f2 = 1.0 - f1

    if last == "1":
        p_add_1 = (r1 * f1) / (r1 * f1 + f2)
    elif last == "2":
        p_add_1 = f1 / (f1 + r2 * f2)
    else:
        p_add_1 = f1

    return "1" if random.random() < p_add_1 else "2"


def simulate_chain(length, f1, r1, r2):
    chain = []

    last = None

    for _ in range(length):
        monomer = choose_next_monomer(last, f1, r1, r2)
        chain.append(monomer)
        last = monomer

    return "".join(chain)


def simulate_chain_with_termination(max_dp, f1, r1, r2, p_terminate):
    chain = []
    last = None

    while len(chain) < max_dp:
        monomer = choose_next_monomer(last, f1, r1, r2)
        chain.append(monomer)
        last = monomer

        if random.random() < p_terminate:
            break

    return "".join(chain)


def block_lengths(chain):
    if not chain:
        return {"1": [], "2": []}

    blocks = {"1": [], "2": []}

    current = chain[0]
    count = 1

    for monomer in chain[1:]:
        if monomer == current:
            count += 1
        else:
            blocks[current].append(count)
            current = monomer
            count = 1

    blocks[current].append(count)
    return blocks


def chain_mass(chain, mw1, mw2):
    return sum(mw1 if unit == "1" else mw2 for unit in chain)


def simulate_copolymerization(
    n_chains=1000,
    target_dp=100,
    f1=0.5,
    r1=1.0,
    r2=1.0,
    mw1=None,
    mw2=None,
    seed=123,
    p_terminate=None,
    max_dp=10000,
):
    random.seed(seed)

    if mw1 is None or mw2 is None:
        raise ValueError("mw1 and mw2 must be provided")

    if p_terminate is None:
        chains = [
            simulate_chain(length=target_dp, f1=f1, r1=r1, r2=r2)
            for _ in range(n_chains)
        ]
    else:
        chains = [
            simulate_chain_with_termination(
                max_dp=max_dp,
                f1=f1,
                r1=r1,
                r2=r2,
                p_terminate=p_terminate,
            )
            for _ in range(n_chains)
        ]

    dps = [len(chain) for chain in chains]
    mean_dp = sum(dps) / len(dps)

    masses = [chain_mass(chain, mw1, mw2) for chain in chains]

    mn = sum(masses) / len(masses)
    mw = sum(m * m for m in masses) / sum(masses)
    
    dispersity = mw / mn

    counts = Counter("".join(chains))
    total_units = counts["1"] + counts["2"]

    fraction_m1 = counts["1"] / total_units
    fraction_m2 = counts["2"] / total_units

    all_blocks_1 = []
    all_blocks_2 = []

    for chain in chains:
        blocks = block_lengths(chain)
        all_blocks_1.extend(blocks["1"])
        all_blocks_2.extend(blocks["2"])

    avg_block_m1 = sum(all_blocks_1) / len(all_blocks_1) if all_blocks_1 else 0.0
    avg_block_m2 = sum(all_blocks_2) / len(all_blocks_2) if all_blocks_2 else 0.0

    return KMCResult(
        chains=chains,
        mn=mn,
        mw=mw,
        dispersity=dispersity,
        mean_dp=mean_dp,
        fraction_m1=fraction_m1,
        fraction_m2=fraction_m2,
        avg_block_m1=avg_block_m1,
        avg_block_m2=avg_block_m2,
    )


if __name__ == "__main__":
    result = simulate_copolymerization(
        n_chains=1000,
        target_dp=100,
        f1=0.5,
        r1=0.395,
        r2=0.586,
        mw1=104.152,
        mw2=130.143,
        seed=123,
        p_terminate=0.01,
        max_dp=10000,
    )

    print(f"Mn: {result.mn:.2f}")
    print(f"Mw: {result.mw:.2f}")
    print(f"Dispersity: {result.dispersity:.3f}")
    print(f"Mean DP: {result.mean_dp:.2f}")
    print(f"Fraction M1: {result.fraction_m1:.3f}")
    print(f"Fraction M2: {result.fraction_m2:.3f}")
    print(f"Average M1 block length: {result.avg_block_m1:.3f}")
    print(f"Average M2 block length: {result.avg_block_m2:.3f}")