project('edipy', 'fortran')


scifor_dep= dependency('scifor', required:true)
edipack_dep= dependency('edipack', required:true)

py = import('python').find_installation(pure: false)

fortran_src = ['src/INPUT_VARS.f90','src/edi2py/edi2py.f90']

library('edi2py', 
                    fortran_src,
                    fortran_args: ['-ffree-line-length-none', '-cpp', '-D_MPI'],
                    dependencies: [scifor_dep,edipack_dep],
                    install: true,
                    install_dir: 'edipy'
)

# Include Python source files
py_files = files('python/edi2py.py', 'python/__init__.py')

install_data(py_files, 
            install_dir:  'edipy'
            )
