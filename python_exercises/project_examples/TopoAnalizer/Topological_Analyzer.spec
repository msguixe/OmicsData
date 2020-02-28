# -*- mode: python -*-

block_cipher = None


a = Analysis(['Topological_Analyzer_to_exe.py'],
             pathex=['C:\\Users\\laure\\AppData\\Local\\Programs\\Python\\Python36-32\\Scripts'],
             binaries=[],
             datas=[],
             hiddenimports=[],
             hookspath=[],
             runtime_hooks=[],
             excludes=[],
             win_no_prefer_redirects=False,
             win_private_assemblies=False,
             cipher=block_cipher)
a.datas += [('help.png','C:\\Users\\laure\\Desktop\\Task\\help.png','DATA'), ('topo_analyzer.ico','C:\\Users\\laure\\Desktop\\Task\\topo_analyzer.ico','DATA') ]
pyz = PYZ(a.pure, a.zipped_data,
             cipher=block_cipher)
exe = EXE(pyz,
          a.scripts,
          a.binaries,
          a.zipfiles,
          a.datas,
          name='Topo_Analyzer',
          debug=False,
          strip=False,
          upx=True,
          runtime_tmpdir=None,
          console=False,
	  icon='C:\\Users\\laure\\Desktop\\topo_analyzer.ico' )
