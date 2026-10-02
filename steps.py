"""Step strings of a study: "A + B <=> C via TS", with optional integer coefficients"""

import re
from dataclasses import dataclass

SEPARATOR = re.compile(r"<=>|⇌|=")
NAME = re.compile(r"[A-Za-z][^\s+*=⇌<>]*")
TERM = re.compile(r"^(?:(?P<coef>\d+)\s*\*?\s+|(?P<coef2>\d+)\*)?(?P<name>.*)$")


class StepSyntaxError(ValueError):
    """A step or overall-reaction string that cannot be parsed"""


@dataclass(frozen=True)
class Step:
    """One reversible elementary step; coefficients are kept per species"""

    reactants: tuple
    products: tuple
    ts: object
    text: str

    def species(self):
        """Reactant and product names, first-seen order, without the transition state"""
        seen = []
        for name, _ in self.reactants + self.products:
            if name not in seen:
                seen.append(name)
        return seen

    def coefficient(self, name, side):
        """Coefficient of name on side ("reactants" or "products"), 0 if absent"""
        return dict(getattr(self, side)).get(name, 0)

    def copasi_equation(self):
        """The step as copasi_helper reads it, coefficients written as repeated species"""

        def side(terms):
            return " + ".join(name for name, n in terms for _ in range(n))

        return f"{side(self.reactants)} = {side(self.products)}"


def _side(text, label, whole):
    if not text.strip():
        raise StepSyntaxError(f'"{whole}": the {label} side is empty')
    counts = {}
    for raw in text.split("+"):
        term = raw.strip()
        if not term:
            raise StepSyntaxError(f'"{whole}": empty species name next to a +')
        match = TERM.match(term)
        coef = match.group("coef") or match.group("coef2")
        name = match.group("name").strip()
        if (
            re.match(r"^\d*\.\d+|^0+\s*\*?\s", term)
            or coef == "0"
            or (coef is not None and int(coef) == 0)
        ):
            raise StepSyntaxError(
                f'"{whole}": coefficients must be positive integers ({term})'
            )
        if " " in name:
            raise StepSyntaxError(f'"{whole}": names cannot contain spaces ({name})')
        if not NAME.fullmatch(name):
            raise StepSyntaxError(f'"{whole}": names must start with a letter ({name})')
        counts[name] = counts.get(name, 0) + (int(coef) if coef else 1)
    return tuple(counts.items())


def parse_step(text, allow_ts=True):
    """Parses "A + 2 B <=> C via TS"; raises StepSyntaxError with the reason"""
    whole = " ".join(str(text).split())
    body, ts = whole, None
    if re.search(r"\svia(\s|$)", whole):
        if not allow_ts:
            raise StepSyntaxError(f'"{whole}": via is not allowed here')
        body, _, after = whole.partition(" via")
        names = after.split()
        if len(names) != 1:
            raise StepSyntaxError(
                f'"{whole}": via must be followed by one transition-state name'
            )
        ts = names[0]
    sides = SEPARATOR.split(body)
    if len(sides) == 1:
        raise StepSyntaxError(f'"{whole}": no <=> between the two sides')
    if len(sides) > 2:
        raise StepSyntaxError(f'"{whole}": more than one <=>')
    left, right = sides
    return Step(_side(left, "left", whole), _side(right, "right", whole), ts, whole)
