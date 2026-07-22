<%namespace module='pyfr.backends.base.makoutil' name='pyfr'/>

<%pyfr:macro name='linear_flux' params='s, f'>
    // f_i = a_i * u   (constant advection velocity dotted with the state)
% for i in range(ndims):
    f[${i}][0] = ${c['a' + 'xyz'[i]]}*s[0];
% endfor
</%pyfr:macro>
