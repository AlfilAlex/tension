# tension

This module allow the user to determine the strechness of a epitope or any other segment of aa. The implementation is simply the use of a distante between aa in a betha barrel, wich usually are streched.
Using information of the Stryer book (Biochemistry), i found the mean distance, so given this distance and the number of aa in a epitope, it is possible to calculate how much a epitope could be streached. Also,, using biipython, is posible to calculate the actual distance between the first an last aa in a epitope chain, so the relative streachnes of a epitope is given by the followinf eequation:

$$S_{rel} = \frac{d_{first,last}}{d_{\beta }*num_{aa}}$$
