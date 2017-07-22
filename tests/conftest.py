import sys
import pkgutil


# def current_module_path(mname):
#     try:
#         return pkgutil.get_loader(mname).filename
#     except AttributeError:
#         return None


# def pytest_cmdline_preparse(config, args):
#     print('root conftest.py: pytest_cmdline_preparse')
#     print('  rif path:', current_module_path('rif'))
#     print('  _rif path:', current_module_path('_rif'))
#     print("  sys.path:")
#     for p in sys.path:
#         print('   ', p)


# def pytest_configure(config):
#     print('root conftest.py: pytest_configure')
#     print('  rif path:', current_module_path('rif'))
#     print('  _rif path:', current_module_path('_rif'))
#     print("  sys.path:")
#     for p in sys.path:
#         print('   ', p)


# def pytest_runtest_setup(item):
#     print ("root conftest.py: setting up", item)
#     print('  rif path:', current_module_path('rif'))
#     print('  _rif path:', current_module_path('_rif'))
#     print("  sys.path:")
#     for p in sys.path:
#         print('   ', p)
