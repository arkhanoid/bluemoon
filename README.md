# bluemoon
Calculate the refraction and absorption of the light through a planetary atmosphere

The atmospheric refraction of an atmosphere is a Luneburg refraction problem, because the index os refraction is spherically simmetric.

To determine the path, whe should solve Euler-Lagrange equations with lagrangian $L=  n(r) * \sqrt{1+(r*d\theta/dr)^2} * dr$

As the light path is calculated, we integrate the absorption coefficient (that is another function of r, due to spheric simetry) over $ds = dr * \sqrt{1+(d\theta/dr)^2}$

Then we apply the Beer lambert law $I = I_0 * e^{- \int a(r) * ds }$ 



