import setuptools

setuptools.setup(
    name="molecule-viewer",
    version="0.0.1",
    packages=setuptools.find_packages(),
    entry_points=dict(console_scripts="molecule-viewer = molecule_viewer.viewer:main")
)