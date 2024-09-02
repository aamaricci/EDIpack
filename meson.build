project('edipy', 'fortran')

python = import('python').find_installation(pure: false)

scifor_dep= dependency('scifor', required:true)
edipack_dep= dependency('edipack', required:true)

fortran_src = ['src/INPUT_VARS.f90','src/edi2py/edi2py.f90']
python_src  = ['python/edi2py.py', 'python/__init__.py']

library('edi2py', 
                    fortran_src,
                    fortran_args: ['-ffree-line-length-none', '-cpp', '-D_MPI'],
                    dependencies: [scifor_dep,edipack_dep],
                    install: true,
                    install_dir: python.get_install_dir() / 'edipy'
)

python.install_sources(
    python_src,
    subdir: 'edipy'
)
