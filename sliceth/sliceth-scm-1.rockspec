package = "sliceth"
version = "scm-1"

source = {
   url = "git@github.com:PhilippPelz/slicepp.git",
   tag = "master"
}

description = {
   summary = "multithreaded multislice code - torch interface",
   detailed = [[
   	    multithreaded multislice code - torch interface
   ]],
   homepage = "https://github.com/PhilippPelz/slicepp"
}

dependencies = {
  "lua >= 5.1",
  "torch >= 7.0",
}

build = {
   type = "command",
   build_command = [[
cmake -E make_directory build && cd build && cmake .. -DCMAKE_BUILD_TYPE=Release -DLUA=$(LUA) -DLUALIB=$(LUALIB) -DLUA_BINDIR="$(LUA_BINDIR)" -DLUA_INCDIR="$(LUA_INCDIR)" -DLUA_LIBDIR="$(LUA_LIBDIR)" -DLUADIR="$(LUADIR)" -DLIBDIR="$(LIBDIR)" -DCMAKE_INSTALL_PREFIX="$(PREFIX)"
   ]],
   install_command = "cd build && make -j$(getconf _NPROCESSORS_ONLN) install"
}
