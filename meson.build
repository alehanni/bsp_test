
project('luajit_test', 'cpp', default_options : ['cpp_std=c++20'])

cpp = meson.get_compiler('cpp')

cmake = import('cmake')
deps = []

# build raylib
rl_opt_var = cmake.subproject_options()
rl_opt_var.add_cmake_defines({
    'BUILD_SHARED_LIBS': 'OFF',
    'USE_EXTERNAL_GLFW': 'OFF'
})

if host_machine.system() == 'emscripten'
    rl_opt_var.add_cmake_defines({'PLATFORM': 'Web'})
endif

raylib_sub_proj = cmake.subproject('raylib', options: rl_opt_var)
raylib_dep = raylib_sub_proj.dependency('raylib')

deps += raylib_dep

if host_machine.system() != 'emscripten'
    deps += raylib_sub_proj.dependency('glfw')
endif

inc = []
sources = []

inc += include_directories('include')
subdir('src')

executable('demo', sources, include_directories : inc, dependencies : deps)
