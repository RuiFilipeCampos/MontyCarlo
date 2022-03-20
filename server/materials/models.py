from django.db import models

# Create your models here.
ELEMENTS = \
        r"H, He, Li, Be, B, C, N, O, F, Ne, Na, Mg, Al, Si, P"      \
        r"S, Cl, Ar, K, Ca, Sc, Ti, V, Cr, Mn, Fe, Co, Ni, Cu, Zn"  \
        r"Ga, Ge, As, Se, Br, Kr, Rb, Sr, Y, Zr, Nb, Mo, Tc, Ru"    \
        r"Rh, Pd, Ag, Cd, In, Sn, Sb, Te, I, Xe, Cs, Ba, La, Ce"    \
        "Pr, Nd, Pm, Sm, Eu, Gd, Tb, Dy, Ho, Er, Tm, Yb, Lu, Hf"   \
        "Ta, W, Re, Os, Ir, Pt, Au, Hg, Tl, Pb, Bi, Po, At, Rn"    \
        "Fr, Ra, Ac, Th, Pa, U, Np, Pu, Am, Cm, Bk, Cf, Es, Fm"    \
        "Md, No, Lr, Rf, Db, Sg, Bh, Hs, Mt, Ds, Rg, Cn, Uut"      \
        "Fl, Uup, Lv, Uus, Uuo"


ELEMENTS = ELEMENTS.split(", ")
CODE = ""
for element_str in ELEMENTS:
    CODE += f"{element_str} = models.IntegerField(blank=True, null=True) \n"


class Material(models.Model):
    name = models.CharField(max_length=500)
    density = models.FloatField()

    state = models.CharField(
        max_length=10,
        blank=True,
        choices=[
            ("LIQUID", "LIQUID"),
            ("SOLID", "SOLID"),
            ("GAS", "GAS"),
        ]
    )

    exec(CODE)








