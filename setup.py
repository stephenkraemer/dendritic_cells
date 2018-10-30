from setuptools import setup, find_packages
setup(
        name='dendritic_cells',
        version='0.1',
        author='Stephen Kraemer',
        author_email='stephenkraemer@gmail.com',
        license='MIT',
        packages=find_packages(),
        # These are the requirements to install the package and introspect it,
        # eg. to get metadata tables.
        # To run reproduce_analysis, or individual scripts or workflows,
        # all packages from the env.yaml must be installed, as described in
        # the reproduce_analysis module
        install_requires=[
            'pandas>=0.23',
            'dss_workflow',
        ],
        python_requires='>=3.6',
        package_data={'': ['*.sh', '*.smk', '*.yml']}
)
