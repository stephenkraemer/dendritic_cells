from setuptools import setup, find_packages
setup(
        name='dendritic_cells',
        version='0.1',
        author='Stephen Kraemer',
        author_email='stephenkraemer@gmail.com',
        license='MIT',
        packages=find_packages(),
        install_requires=[
            'pandas>=0.23',
            'snakemake>=5.3.0',
        ],
        python_requires='>=3.6',
        package_data={'': ['*.sh', '*.smk', '*.yml']}
)
