# FastSpie 4

The FastSpie4 is similar to the FastSpie3 in the splitting method, but uses a better error estimation, yielding complexity improvements.

The previous FastSpies computed error bounds using the worst possible estimate: The positive error bound is that all products were created immediately, and the negative is that all products were eliminated immediately.

Instead, we can sample the actual timings of the events. For $k$ work, we get a precision of $2^k$, which is not bad.

What does this mean for our bounds?

When a reactant is rising, it does nothing: The old error bound is not accurate, but it is half accurate. However, when a reaction is generally constant, it does reduce the errors to a square root of what they used to be.

We can do better.

To fix this, we will also estimate the derivative, and sample the number of events using the first derivative. We are left with two errors, which I will show on the $C \to D$ part of $A \stackrel{r_1}{\to} B \stackrel{r_2}{\to} C \stackrel{r_3}{\to} D$:

* The third-order error from second-derivative effects. For the above reaction these are $O([A]r_1r_2r_3t^3)$.
* The errors from the stochastic deviation. These are $O\left(\sqrt{[B]r_2t}r_3t\right) = O\left([B]^{0.5}r_2^{0.5}r_3t^{1.5}\right)$.

Normalizing both errors to be 1 and the number of reactants to be $N$, we get $t \sim N^{-\frac{1}{3}}$.

Also note that the stochastic deviation errors are distributed normally, so they sum neatly when multiple reactions are composed. The third-order error doesn't, but it is an inherent property of the SDE and isn't affected by a coordinate change.

For $N=2^{15}$ as in the BCR system I am not sure the difference is that large. The reaction averaging effect would be nice.