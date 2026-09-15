# -*- mode: python ; coding: utf-8 -*-
import os

# Add src to the module search path
pathex = [os.path.abspath('src')]

# Bundle static assets (no more worrying about Windows/Mac path-separator differences)
datas = [
    ('src/castorGUI/frontend', 'frontend'),
    ('src/castorGUI/data', 'data'),
]

# List every module to exclude here; add as many as you like.
#
# astropy.visualization pulls in matplotlib at import time (wcsaxes/__init__.py
# calls pytest.importorskip("matplotlib") itself), and CASTOR only ever touches
# astropy.time / astropy.coordinates / astropy.units — so it goes too, along
# with everything matplotlib alone was dragging in (Pillow, fonttools, ...).
excludes = [
    'pytest',
    'matplotlib',
    'astropy.visualization',
    'tkinter',
    'IPython',
    'notebook'
]

a = Analysis(
    ['src/castorGUI/desktop.py'],
    pathex=pathex,
    binaries=[],
    datas=datas,
    hiddenimports=[],
    # Shadows pyinstaller-hooks-contrib's hook-astropy.py, which crashes the
    # whole build over one submodule needing matplotlib — see the docstring
    # in pyinstaller_hooks/hook-astropy.py.
    hookspath=['pyinstaller_hooks'],
    hooksconfig={},
    runtime_hooks=[],
    excludes=excludes,
    noarchive=False,
)

pyz = PYZ(a.pure)

exe = EXE(
    pyz,
    a.scripts,
    a.binaries,
    a.datas,
    [],
    name='CASTOR-ETC', # name of the output executable
    icon='assets/desktop/castor.ico', # Windows executable icon
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=True,
    upx_exclude=[],
    runtime_tmpdir=None,
    console=False, # equivalent to the original --noconsole
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
)

app = BUNDLE(
    exe,
    name='CASTOR-ETC.app',
    icon='assets/desktop/castor.icns',
    bundle_identifier=None,
)
