function out = model
%
% aniso.m
%
% Model exported on Feb 5 2019, 15:03 by COMSOL 5.3.1.275.

import com.comsol.model.*
import com.comsol.model.util.*

model = ModelUtil.create('Model');

model.modelPath('/home/vd256574/Formation_TRUST/poisson_anisotrope/comsol');

model.label('aniso.mph');

model.param.set('x1', '0.5');
model.param.set('x2', '1.5');
model.param.set('a', '100', 'diffusivite horizontale');
model.param.set('b', '1', 'diffusivite verticale');
model.param.set('stepsize', '0.2');

model.component.create('comp1', true);

model.component('comp1').geom.create('geom1', 2);

model.component('comp1').label('Composant 1');

model.result.table.create('tbl1', 'Table');
model.result.table.create('tbl2', 'Table');
model.result.table.create('tbl3', 'Table');

model.component('comp1').func.create('pw1', 'Piecewise');
model.component('comp1').func.create('pw2', 'Piecewise');
model.component('comp1').func('pw1').label('SmoothStep');
model.component('comp1').func('pw1').set('funcname', 'stepf');
model.component('comp1').func('pw1').set('pieces', {'x1-0.5*stepsize' 'x1+0.5*stepsize' '-sin((x-x1)*pi/stepsize)'; 'x1+stepsize*0.5' 'x2-stepsize*0.5' '-1'; 'x2-0.5*stepsize' 'x2+0.5*stepsize' 'sin((x-x2)*pi/stepsize)'});
model.component('comp1').func('pw2').set('pieces', {'0' '0.4' '1'; '0.4' '1' '0'});

model.component('comp1').mesh.create('mesh1');

model.component('comp1').geom('geom1').label(['G' native2unicode(hex2dec({'00' 'e9'}), 'unicode') 'om' native2unicode(hex2dec({'00' 'e9'}), 'unicode') 'trie 1']);
model.component('comp1').geom('geom1').create('pc1', 'ParametricCurve');
model.component('comp1').geom('geom1').feature('pc1').label(['Courbe param' native2unicode(hex2dec({'00' 'e9'}), 'unicode') 'tr' native2unicode(hex2dec({'00' 'e9'}), 'unicode') 'e 1']);
model.component('comp1').geom('geom1').feature('pc1').set('coord', {'10*s' 'stepf(s)*0.25+0.75'});
model.component('comp1').geom('geom1').create('ls1', 'LineSegment');
model.component('comp1').geom('geom1').feature('ls1').label(['Segment lin' native2unicode(hex2dec({'00' 'e9'}), 'unicode') 'aire 1']);
model.component('comp1').geom('geom1').feature('ls1').set('specify1', 'coord');
model.component('comp1').geom('geom1').feature('ls1').set('coord1', [0 1]);
model.component('comp1').geom('geom1').feature('ls1').set('specify2', 'coord');
model.component('comp1').geom('geom1').create('ls2', 'LineSegment');
model.component('comp1').geom('geom1').feature('ls2').label(['Segment lin' native2unicode(hex2dec({'00' 'e9'}), 'unicode') 'aire 2']);
model.component('comp1').geom('geom1').feature('ls2').set('specify1', 'coord');
model.component('comp1').geom('geom1').feature('ls2').set('specify2', 'coord');
model.component('comp1').geom('geom1').feature('ls2').set('coord2', [10 0]);
model.component('comp1').geom('geom1').create('ls3', 'LineSegment');
model.component('comp1').geom('geom1').feature('ls3').label(['Segment lin' native2unicode(hex2dec({'00' 'e9'}), 'unicode') 'aire 3']);
model.component('comp1').geom('geom1').feature('ls3').set('specify1', 'coord');
model.component('comp1').geom('geom1').feature('ls3').set('coord1', [10 0]);
model.component('comp1').geom('geom1').feature('ls3').set('specify2', 'coord');
model.component('comp1').geom('geom1').feature('ls3').set('coord2', [10 1]);
model.component('comp1').geom('geom1').create('csol1', 'ConvertToSolid');
model.component('comp1').geom('geom1').feature('csol1').label('Convertir en solide 1');
model.component('comp1').geom('geom1').feature('csol1').selection('input').set({'ls1' 'ls2' 'ls3' 'pc1'});
model.component('comp1').geom('geom1').feature('fin').label('Constituer une union');
model.component('comp1').geom('geom1').run;
model.component('comp1').geom('geom1').run('fin');

model.component('comp1').variable.create('var1');
model.component('comp1').variable('var1').set('c', 'u2y/sqrt(u2x*u2x+u2y*u2y)');
model.component('comp1').variable('var1').set('s', 'u2x/sqrt(u2x*u2x+u2y*u2y)');
model.component('comp1').variable('var1').set('d11', 'a*c*c+b*s*s', 'D[1][1]');
model.component('comp1').variable('var1').set('d12', 'c*s*(b-a)', 'D[1][2] ou D[2][1]');
model.component('comp1').variable('var1').set('d22', 'a*s*s+b*c*c', 'D[2][2]');
model.component('comp1').variable.create('var2');
model.component('comp1').variable('var2').set('q11', '-(d11*ux+d12*uy)');
model.component('comp1').variable('var2').set('q12', '-(d12*ux+d22*uy)');
model.component('comp1').variable('var2').set('q', 'sqrt(q11*q11+q12*q12)');
model.component('comp1').variable('var2').selection.geom('geom1', 1);
model.component('comp1').variable('var2').selection.set([4]);

model.view.create('view2', 2);
model.view.create('view3', 2);

model.component('comp1').physics.create('lpeq', 'LaplaceEquation', 'geom1');
model.component('comp1').physics('lpeq').field('dimensionless').field('u2');
model.component('comp1').physics('lpeq').create('dir1', 'DirichletBoundary', 1);
model.component('comp1').physics('lpeq').feature('dir1').selection.set([2]);
model.component('comp1').physics('lpeq').create('dir2', 'DirichletBoundary', 1);
model.component('comp1').physics('lpeq').feature('dir2').selection.set([4]);
model.component('comp1').physics.create('poeq', 'PoissonEquation', 'geom1');
model.component('comp1').physics('poeq').create('flux1', 'FluxBoundary', 1);
model.component('comp1').physics('poeq').feature('flux1').selection.set([4]);
model.component('comp1').physics('poeq').create('dir1', 'DirichletBoundary', 1);
model.component('comp1').physics('poeq').feature('dir1').selection.set([2]);
model.component('comp1').physics('poeq').create('cons1', 'Constraint', 1);
model.component('comp1').physics('poeq').feature('cons1').selection.set([4]);
model.component('comp1').physics('poeq').create('disc1', 'Discretization', -1);

model.component('comp1').mesh('mesh1').create('ftri1', 'FreeTri');

model.result.table('tbl1').comments(['Int' native2unicode(hex2dec({'00' 'e9'}), 'unicode') 'gration sur ligne 1 (dflux.u, ux*nx+uy*ny)']);
model.result.table('tbl2').comments(['Int' native2unicode(hex2dec({'00' 'e9'}), 'unicode') 'gration sur ligne 1 (dflux.u, ux*nx+uy*ny)']);
model.result.table('tbl3').comments(['Int' native2unicode(hex2dec({'00' 'e9'}), 'unicode') 'gration sur ligne 1 (dflux.u, ux*nx+uy*ny)']);

model.capeopen.label('Thermodynamique');

model.component('comp1').variable('var2').label('flux_sur_surface');

model.component('comp1').view('view1').label('Vue 1');
model.component('comp1').view('view1').set('showlabels', true);
model.component('comp1').view('view1').set('showDirections', true);
model.component('comp1').view('view1').axis.label('Axe');
model.component('comp1').view('view1').axis.set('xmin', -1.901564121246338);
model.component('comp1').view('view1').axis.set('xmax', 11.90156364440918);
model.component('comp1').view('view1').axis.set('ymin', -3.9880404472351074);
model.component('comp1').view('view1').axis.set('ymax', 4.988040447235107);
model.component('comp1').view('view1').axis.set('abstractviewlratio', -0.05000000074505806);
model.component('comp1').view('view1').axis.set('abstractviewrratio', 0.05000000074505806);
model.component('comp1').view('view1').axis.set('abstractviewbratio', -3.79068660736084);
model.component('comp1').view('view1').axis.set('abstractviewtratio', 3.7906863689422607);
model.component('comp1').view('view1').axis.set('abstractviewxscale', 0.010119595564901829);
model.component('comp1').view('view1').axis.set('abstractviewyscale', 0.010119596496224403);
model.view('view2').label('Vue 2D 2');
model.view('view2').axis.label('Axe');
model.view('view2').axis.set('xmin', -6.000001430511475);
model.view('view2').axis.set('xmax', 16);
model.view('view2').axis.set('ymin', -0.6000032424926758);
model.view('view2').axis.set('ymax', 1.600008487701416);
model.view('view2').axis.set('viewscaletype', 'automatic');
model.view('view2').axis.set('abstractviewlratio', -0.6000001430511475);
model.view('view2').axis.set('abstractviewrratio', 0.6000000238418579);
model.view('view2').axis.set('abstractviewbratio', -0.6000000834465027);
model.view('view2').axis.set('abstractviewtratio', 0.6000000834465027);
model.view('view2').axis.set('abstractviewxscale', 0.020618557929992676);
model.view('view2').axis.set('abstractviewyscale', 0.0025974165182560682);
model.view('view3').axis.set('xmin', -6.000001430511475);
model.view('view3').axis.set('xmax', 16);
model.view('view3').axis.set('ymin', -0.6000032424926758);
model.view('view3').axis.set('ymax', 1.600008487701416);
model.view('view3').axis.set('viewscaletype', 'automatic');
model.view('view3').axis.set('abstractviewlratio', -0.6000001430511475);
model.view('view3').axis.set('abstractviewrratio', 0.6000000238418579);
model.view('view3').axis.set('abstractviewbratio', -0.6000000834465027);
model.view('view3').axis.set('abstractviewtratio', 0.6000000834465027);
model.view('view3').axis.set('abstractviewxscale', 0.020618557929992676);
model.view('view3').axis.set('abstractviewyscale', 0.0025974165182560682);

model.component('comp1').coordSystem('sys1').label(['Rep' native2unicode(hex2dec({'00' 'e8'}), 'unicode') 're sur fronti' native2unicode(hex2dec({'00' 'e8'}), 'unicode') 're 1']);

model.component('comp1').physics('lpeq').label('Equation de Laplace');
model.component('comp1').physics('lpeq').feature('leq1').label('Equation de Laplace 1');
model.component('comp1').physics('lpeq').feature('leq1').featureInfo('info').label(['Vue des ' native2unicode(hex2dec({'00' 'e9'}), 'unicode') 'quations']);
model.component('comp1').physics('lpeq').feature('zflx1').label('Flux nul 1');
model.component('comp1').physics('lpeq').feature('zflx1').featureInfo('info').label(['Vue des ' native2unicode(hex2dec({'00' 'e9'}), 'unicode') 'quations']);
model.component('comp1').physics('lpeq').feature('init1').set('u2', 'y');
model.component('comp1').physics('lpeq').feature('init1').label('Valeurs initiales 1');
model.component('comp1').physics('lpeq').feature('init1').featureInfo('info').label(['Vue des ' native2unicode(hex2dec({'00' 'e9'}), 'unicode') 'quations']);
model.component('comp1').physics('lpeq').feature('dir1').label('Condition de Dirichlet 1');
model.component('comp1').physics('lpeq').feature('dir1').featureInfo('info').label(['Vue des ' native2unicode(hex2dec({'00' 'e9'}), 'unicode') 'quations']);
model.component('comp1').physics('lpeq').feature('dir2').set('r', 1);
model.component('comp1').physics('lpeq').feature('dir2').label('Condition de Dirichlet 2');
model.component('comp1').physics('lpeq').feature('dir2').featureInfo('info').label(['Vue des ' native2unicode(hex2dec({'00' 'e9'}), 'unicode') 'quations']);
model.component('comp1').physics('poeq').label('Equation de Poisson');
model.component('comp1').physics('poeq').feature('peq1').set('f', 0);
model.component('comp1').physics('poeq').feature('peq1').set('c', {'a*c*c+b*s*s' 'c*s*(b-a)' 'c*s*(b-a)' 'a*s*s+b*c*c'});
model.component('comp1').physics('poeq').feature('peq1').label('Equation de Poisson 1');
model.component('comp1').physics('poeq').feature('peq1').featureInfo('info').label(['Vue des ' native2unicode(hex2dec({'00' 'e9'}), 'unicode') 'quations']);
model.component('comp1').physics('poeq').feature('zflx1').label('Flux nul 1');
model.component('comp1').physics('poeq').feature('zflx1').featureInfo('info').label(['Vue des ' native2unicode(hex2dec({'00' 'e9'}), 'unicode') 'quations']);
model.component('comp1').physics('poeq').feature('init1').label('Valeurs initiales 1');
model.component('comp1').physics('poeq').feature('init1').featureInfo('info').label(['Vue des ' native2unicode(hex2dec({'00' 'e9'}), 'unicode') 'quations']);
model.component('comp1').physics('poeq').feature('flux1').set('g', '-0.5*(1-stepf(x/10))*(u-1)*10');
model.component('comp1').physics('poeq').feature('flux1').label('Flux ou source 1');
model.component('comp1').physics('poeq').feature('flux1').featureInfo('info').label(['Vue des ' native2unicode(hex2dec({'00' 'e9'}), 'unicode') 'quations']);
model.component('comp1').physics('poeq').feature('dir1').label('Condition de Dirichlet 1');
model.component('comp1').physics('poeq').feature('dir1').featureInfo('info').label(['Vue des ' native2unicode(hex2dec({'00' 'e9'}), 'unicode') 'quations']);
model.component('comp1').physics('poeq').feature('cons1').set('R', '(1-stepf(x/10))*10000*(1-u)-q*0.5*(1+stepf(x/10))');
model.component('comp1').physics('poeq').feature('cons1').active(false);
model.component('comp1').physics('poeq').feature('cons1').label('Contrainte 1');
model.component('comp1').physics('poeq').feature('cons1').featureInfo('info').label(['Vue des ' native2unicode(hex2dec({'00' 'e9'}), 'unicode') 'quations']);
model.component('comp1').physics('poeq').feature('disc1').label(['Discr' native2unicode(hex2dec({'00' 'e9'}), 'unicode') 'tisation 1']);
model.component('comp1').physics('poeq').feature('disc1').featureInfo('info').label(['Vue des ' native2unicode(hex2dec({'00' 'e9'}), 'unicode') 'quations']);

model.component('comp1').mesh('mesh1').label('Maillage 1');
model.component('comp1').mesh('mesh1').feature('size').label('Taille');
model.component('comp1').mesh('mesh1').feature('size').set('hauto', 1);
model.component('comp1').mesh('mesh1').feature('size').set('custom', 'on');
model.component('comp1').mesh('mesh1').feature('size').set('hmax', 0.05);
model.component('comp1').mesh('mesh1').run;

model.study.create('std1');
model.study('std1').create('stat', 'Stationary');

model.sol.create('sol1');
model.sol('sol1').study('std1');
model.sol('sol1').attach('std1');
model.sol('sol1').create('st1', 'StudyStep');
model.sol('sol1').create('v1', 'Variables');
model.sol('sol1').create('s1', 'Stationary');
model.sol('sol1').feature('s1').create('fc1', 'FullyCoupled');
model.sol('sol1').feature('s1').feature.remove('fcDef');

model.result.dataset.create('edg1', 'Edge2D');
model.result.dataset.create('cln1', 'CutLine2D');
model.result.dataset('edg1').selection.set([4]);
model.result.numerical.create('int1', 'IntLine');
model.result.numerical.create('int2', 'IntLine');
model.result.numerical('int1').selection.set([4]);
model.result.numerical('int1').set('probetag', 'none');
model.result.numerical('int2').selection.set([2]);
model.result.numerical('int2').set('probetag', 'none');
model.result.create('pg1', 'PlotGroup2D');
model.result.create('pg2', 'PlotGroup1D');
model.result.create('pg4', 'PlotGroup1D');
model.result.create('pg5', 'PlotGroup2D');
model.result('pg1').create('surf1', 'Surface');
model.result('pg1').create('con1', 'Contour');
model.result('pg2').create('lngr1', 'LineGraph');
model.result('pg2').create('lngr2', 'LineGraph');
model.result('pg2').create('lngr3', 'LineGraph');
model.result('pg2').create('lngr4', 'LineGraph');
model.result('pg2').create('lngr5', 'LineGraph');
model.result('pg2').feature('lngr1').set('xdata', 'expr');
model.result('pg2').feature('lngr2').set('xdata', 'expr');
model.result('pg2').feature('lngr2').selection.set([2]);
model.result('pg2').feature('lngr3').set('xdata', 'expr');
model.result('pg4').create('lngr1', 'LineGraph');
model.result('pg5').create('surf1', 'Surface');
model.result('pg5').create('surf2', 'Surface');
model.result('pg5').create('surf3', 'Surface');
model.result.export.create('plot2', 'Plot');
model.result.export.create('plot3', 'Plot');

model.study('std1').label('Etude 1');
model.study('std1').feature('stat').label('Stationnaire');

model.sol('sol1').attach('std1');
model.sol('sol1').runAll;

model.result.label(['R' native2unicode(hex2dec({'00' 'e9'}), 'unicode') 'sultats']);
model.result.dataset('edg1').label('Bord Haut');
model.result.dataset('cln1').label('Diagonal');
model.result.dataset('cln1').set('genpoints', [0 1; 10 0]);
model.result.numerical('int1').label(['Int' native2unicode(hex2dec({'00' 'e9'}), 'unicode') 'gration sur ligne 1']);
model.result.numerical('int1').set('table', 'tbl3');
model.result.numerical('int1').set('expr', {'dflux.u' 'ux*nx+uy*ny' ''});
model.result.numerical('int1').set('unit', {'1' '1' ''});
model.result.numerical('int1').set('descr', {['Flux sur fronti' native2unicode(hex2dec({'00' 'e8'}), 'unicode') 're, direction down'] '' ''});
model.result.numerical('int1').set('dataseries', 'integral');
model.result.numerical('int2').label(['Int' native2unicode(hex2dec({'00' 'e9'}), 'unicode') 'gration sur ligne 2']);
model.result.numerical('int2').set('table', 'tbl3');
model.result.numerical('int2').set('expr', {'uy'});
model.result.numerical('int2').set('unit', {'1'});
model.result.numerical('int2').set('descr', {'Gradient de u, composante y'});
model.result.numerical('int1').setResult;
model.result.numerical('int2').appendResult;
model.result('pg1').label('sonde_dom');
model.result('pg1').feature('surf1').label('temperature');
model.result('pg1').feature('surf1').set('expr', 'u');
model.result('pg1').feature('surf1').set('descr', 'Dependent variable u');
model.result('pg1').feature('surf1').set('resolution', 'normal');
model.result('pg1').feature('con1').label('Isovaleurs temp');
model.result('pg1').feature('con1').set('expr', 'u');
model.result('pg1').feature('con1').set('descr', ['Variable d' native2unicode(hex2dec({'00' 'e9'}), 'unicode') 'pendante u']);
model.result('pg1').feature('con1').set('number', 10);
model.result('pg1').feature('con1').set('coloring', 'uniform');
model.result('pg1').feature('con1').set('color', 'black');
model.result('pg1').feature('con1').set('resolution', 'finer');
model.result('pg1').feature('con1').set('smooth', 'none');
model.result('pg1').feature('con1').set('recover', 'pprint');
model.result('pg1').feature('con1').set('resolution', 'finer');
model.result('pg2').label('sonde_Bord_Haut');
model.result('pg2').set('data', 'edg1');
model.result('pg2').set('xlabel', 'x');
model.result('pg2').set('xlabelactive', true);
model.result('pg2').set('ylabel', ['Flux sur fronti' native2unicode(hex2dec({'00' 'e8'}), 'unicode') 're, direction down (1/m)']);
model.result('pg2').set('ylabelactive', false);
model.result('pg2').feature('lngr1').active(false);
model.result('pg2').feature('lngr1').label('uy*ny+ux*nx');
model.result('pg2').feature('lngr1').set('data', 'edg1');
model.result('pg2').feature('lngr1').set('expr', 'uy*ny+ux*nx');
model.result('pg2').feature('lngr1').set('unit', '1/m');
model.result('pg2').feature('lngr1').set('descr', 'uy*ny+ux*nx');
model.result('pg2').feature('lngr1').set('xdataexpr', 'x');
model.result('pg2').feature('lngr1').set('xdataunit', 'm');
model.result('pg2').feature('lngr1').set('xdatadescr', ['Coordonn' native2unicode(hex2dec({'00' 'e9'}), 'unicode') 'e x']);
model.result('pg2').feature('lngr1').set('linecolor', 'blue');
model.result('pg2').feature('lngr1').set('linewidth', 2);
model.result('pg2').feature('lngr1').set('legend', true);
model.result('pg2').feature('lngr1').set('legendmethod', 'manual');
model.result('pg2').feature('lngr1').set('legends', {'uy*ny+ux*nx'});
model.result('pg2').feature('lngr1').set('resolution', 'normal');
model.result('pg2').feature('lngr2').active(false);
model.result('pg2').feature('lngr2').label('condlim');
model.result('pg2').feature('lngr2').set('data', 'edg1');
model.result('pg2').feature('lngr2').set('expr', '-0.5*(1-stepf(x/10))*(u-1)*10000');
model.result('pg2').feature('lngr2').set('unit', '');
model.result('pg2').feature('lngr2').set('descr', '-0.5*(1-stepf(x/10))*(u-1)*10000');
model.result('pg2').feature('lngr2').set('xdataexpr', 'x');
model.result('pg2').feature('lngr2').set('xdataunit', 'm');
model.result('pg2').feature('lngr2').set('xdatadescr', ['Coordonn' native2unicode(hex2dec({'00' 'e9'}), 'unicode') 'e x']);
model.result('pg2').feature('lngr2').set('linecolor', 'red');
model.result('pg2').feature('lngr2').set('linewidth', 2);
model.result('pg2').feature('lngr2').set('legend', true);
model.result('pg2').feature('lngr2').set('legendmethod', 'manual');
model.result('pg2').feature('lngr2').set('legends', {'uy'});
model.result('pg2').feature('lngr2').set('resolution', 'normal');
model.result('pg2').feature('lngr3').label('dflux.u');
model.result('pg2').feature('lngr3').set('data', 'edg1');
model.result('pg2').feature('lngr3').set('expr', 'dflux.u');
model.result('pg2').feature('lngr3').set('unit', '1/m');
model.result('pg2').feature('lngr3').set('descr', ['Flux sur fronti' native2unicode(hex2dec({'00' 'e8'}), 'unicode') 're, direction down']);
model.result('pg2').feature('lngr3').set('xdataexpr', 'x');
model.result('pg2').feature('lngr3').set('xdataunit', 'm');
model.result('pg2').feature('lngr3').set('xdatadescr', ['Coordonn' native2unicode(hex2dec({'00' 'e9'}), 'unicode') 'e x']);
model.result('pg2').feature('lngr3').set('linecolor', 'magenta');
model.result('pg2').feature('lngr3').set('linewidth', 2);
model.result('pg2').feature('lngr3').set('legend', true);
model.result('pg2').feature('lngr3').set('legendmethod', 'manual');
model.result('pg2').feature('lngr3').set('legends', {'dflux.u'});
model.result('pg2').feature('lngr3').set('resolution', 'normal');
model.result('pg2').feature('lngr4').active(false);
model.result('pg2').feature('lngr4').label('q');
model.result('pg2').feature('lngr4').set('data', 'edg1');
model.result('pg2').feature('lngr4').set('expr', 'q');
model.result('pg2').feature('lngr4').set('unit', '1/m');
model.result('pg2').feature('lngr4').set('descr', '');
model.result('pg2').feature('lngr4').set('linestyle', 'dashdot');
model.result('pg2').feature('lngr4').set('linecolor', 'black');
model.result('pg2').feature('lngr4').set('linewidth', 2);
model.result('pg2').feature('lngr4').set('legend', true);
model.result('pg2').feature('lngr4').set('legendmethod', 'manual');
model.result('pg2').feature('lngr4').set('legends', {'-q'});
model.result('pg2').feature('lngr4').set('resolution', 'normal');
model.result('pg2').feature('lngr5').active(false);
model.result('pg2').feature('lngr5').label('stepf');
model.result('pg2').feature('lngr5').set('expr', '(1-stepf(x/10))*0.5');
model.result('pg2').feature('lngr5').set('unit', '');
model.result('pg2').feature('lngr5').set('descr', '(1-stepf(x/10))*0.5');
model.result('pg2').feature('lngr5').set('resolution', 'normal');
model.result('pg4').label('sonde_diagonal');
model.result('pg4').set('data', 'cln1');
model.result('pg4').set('xlabel', 'x');
model.result('pg4').set('xlabelactive', true);
model.result('pg4').set('ylabel', 'Dependent variable u (1)');
model.result('pg4').set('ylabelactive', false);
model.result('pg4').feature('lngr1').label('temperature');
model.result('pg4').feature('lngr1').set('expr', 'u');
model.result('pg4').feature('lngr1').set('descr', 'Dependent variable u');
model.result('pg4').feature('lngr1').set('resolution', 'normal');
model.result('pg5').label('coefficient_anisotrope');
model.result('pg5').set('window', 'window4');
model.result('pg5').set('windowtitle', 'Plot 4');
model.result('pg5').feature('surf1').label('D11');
model.result('pg5').feature('surf1').set('expr', 'd11');
model.result('pg5').feature('surf1').set('descr', 'D[1][1]');
model.result('pg5').feature('surf1').set('resolution', 'normal');
model.result('pg5').feature('surf2').label('D12&D21');
model.result('pg5').feature('surf2').set('expr', 'd12');
model.result('pg5').feature('surf2').set('descr', 'D[1][2] ou D[2][1]');
model.result('pg5').feature('surf2').set('resolution', 'normal');
model.result('pg5').feature('surf3').label('D22');
model.result('pg5').feature('surf3').set('expr', 'd22');
model.result('pg5').feature('surf3').set('descr', 'D[2][2]');
model.result('pg5').feature('surf3').set('resolution', 'normal');
model.result.export('plot2').label('Plot_flux_bord_haut');
model.result.export('plot2').set('plotgroup', 'pg2');
model.result.export('plot2').set('plot', 'lngr3');
model.result.export('plot2').set('filename', 'comsol_flux_sur_surface_bord_haut.txt');
model.result.export('plot2').set('header', false);
model.result.export('plot2').set('fullprec', false);
model.result.export('plot3').label('Plot_temperature_diagonal');
model.result.export('plot3').set('plotgroup', 'pg2');
model.result.export('plot3').set('plot', 'lngr3');
model.result.export('plot3').set('filename', 'comsol_temperature_sur_diagonal.txt');
model.result.export('plot3').set('header', false);
model.result.export('plot3').set('fullprec', false);

out = model;
