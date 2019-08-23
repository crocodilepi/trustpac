# First test of implicit time scheme
Fully implicit laplacian with mixed boundary conditions:
$$\int_{E_i}-\nabla\cdot k\nabla x^{t+1} = rhs$$
$rhs$ contains contributions from the source term and mixed boundary condition.
For element $i$ on the boundary, this is added to matrix and rhs for the Boundary
face:
$$a_{ii}x_i^{t+1} + ... = a_{ii}x_{ext}$$
where $x_{ext}$ is the boundary value imposed for x.

The explicit operator computes
$$explicit_operator(x) = a_{ii}(x_i^{t} + ... - a_{ii}x_{ext}$$
Hence the $rhs$ that we need can be obtained by computing
$rhs = explicit_operator(x) + M x^t$  where $M$ is the operator matrix.

This is the what equation::assembler does:
'''
void  Equation_base::assembler( Matrice_Morse& matrice,const DoubleTab& inco, DoubleTab& resu)
{
  for (int op=0; op< nombre_d_operateurs(); op++)
    {
      operateur(op).l_op_base().contribuer_a_avec(inco, matrice );
      operateur(op).ajouter(resu);
    }
  sources().ajouter(resu);
  sources().contribuer_a_avec(inco,matrice);
  matrice.ajouter_multvect(inco,resu);
}
'''
