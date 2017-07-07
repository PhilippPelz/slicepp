<div class="section" id="welcome-to-slicepp-s-documentation">
<h1>Welcome to Slicepp&#8217;s documentation!<a class="headerlink" href="#welcome-to-slicepp-s-documentation" title="Permalink to this headline">¶</a></h1>
<div class="toctree-wrapper compound">
<ul class="simple">
</ul>
</div>
</div>
<div class="section" id="about">
<h1>About<a class="headerlink" href="#about" title="Permalink to this headline">¶</a></h1>
<p>Slicepp is a multislice simulation software.It simulates an incoming electron wave striking and propagating through a user-defined sample, eventually picked up by a detector. This software is accelerated by the NVIDIA CUDA platform and ArrayFire.</p>
</div>
<div class="section" id="installation">
<h1>Installation<a class="headerlink" href="#installation" title="Permalink to this headline">¶</a></h1>
<p><strong>Libraries Required</strong></p>
<table border="1" class="docutils">
<colgroup>
<col width="65%" />
<col width="18%" />
<col width="17%" />
</colgroup>
<thead valign="bottom">
<tr class="row-odd"><th class="head">Library</th>
<th class="head">Version</th>
<th class="head">Comments</th>
</tr>
</thead>
<tbody valign="top">
<tr class="row-even"><td><a class="reference external" href="http://arma.sourceforge.net/">Armadillo</a></td>
<td>&gt;= 5.0</td>
<td>Required</td>
</tr>
<tr class="row-odd"><td><a class="reference external" href="http://www.arrayfire.com/">ArrayFire</a></td>
<td>&gt;= 3.0</td>
<td>Required</td>
</tr>
<tr class="row-even"><td><a class="reference external" href="http://www.boost.org/">Boost</a></td>
<td>&gt;= 1.5</td>
<td>Required</td>
</tr>
<tr class="row-odd"><td><a class="reference external" href="https://developer.nvidia.com/cuda-downloads">CUDA Toolkit</a></td>
<td>&gt;= 7.0</td>
<td>Required</td>
</tr>
<tr class="row-even"><td><a class="reference external" href="http://www.fftw.org/">FFTW3</a></td>
<td>&gt;= 3.3</td>
<td>Required</td>
</tr>
<tr class="row-odd"><td><a class="reference external" href="https://www.hdfgroup.org/HDF5/">HDF5</a></td>
<td>&gt;= 1.8</td>
<td>Required</td>
</tr>
<tr class="row-even"><td><a class="reference external" href="http://ab-initio.mit.edu/wiki/index.php/NLopt">NLopt</a></td>
<td>&gt;= 2.4</td>
<td>Required</td>
</tr>
<tr class="row-odd"><td><a class="reference external" href="http://openbabel.org/wiki/Main_Page">OpenBabel</a></td>
<td>&gt;= 2.3</td>
<td>Required</td>
</tr>
<tr class="row-even"><td><a class="reference external" href="http://openmp.org/wp/">OpenMP</a></td>
<td>&gt;= 4.0</td>
<td>Required</td>
</tr>
<tr class="row-odd"><td><a class="reference external" href="https://www.python.org/">Python2</a></td>
<td>2.7</td>
<td>GUI</td>
</tr>
</tbody>
</table>
<p><strong>Compilation</strong></p>
<p><em>Debug</em></p>
<div class="highlight-python"><div class="highlight"><pre>% cd /slicepp/build/debug
% ./remake.sh
% make -j2
</pre></div>
</div>
<p><em>Release</em></p>
<div class="highlight-python"><div class="highlight"><pre>% cd /slicepp/build/release
% ./remake.sh:
</pre></div>
</div>
</div>
<div class="section" id="execution">
<h1>Execution<a class="headerlink" href="#execution" title="Permalink to this headline">¶</a></h1>
<p><em>From terminal</em></p>
<p>Run the following commands in terminal:</p>
<div class="highlight-python"><div class="highlight"><pre>% cd /slicepp/build/release/bin
%./stem3 [PATH_TO_CONFIG_FILE]
</pre></div>
</div>
<p><em>From GUI</em></p>
<p>Run the following commands in terminal:</p>
<div class="highlight-python"><div class="highlight"><pre>% cd /slicepp/GUI_new
% python window.py
</pre></div>
</div>
</div>
<div class="section" id="programming-guide">
<h1>Programming Guide<a class="headerlink" href="#programming-guide" title="Permalink to this headline">¶</a></h1>
<p><strong>Source Code Structure</strong></p>
<p>Source code for the main project are in the <strong>/libs</strong> folder. Families of hierarchial code are grouped into subdirectories (e.g. wavefunctions, potentials, detectors). In each of these subdirectories, there is a header file, usually called <strong>[DIRECTORYNAME]_interface.hpp</strong>, that serves as a template for other classes in the same directory.To create a new class, simply inherit from the template and override methods where appropriate. For instance, if one wants to create a new way of calculate potential, one can make a subclass of potential_interface called NewPotential. Then, implement pure virtual and other additional methods in the NewPotential class. When creating a new subdirectory of <strong>/libs</strong>, make sure also to create a CMakeList.txt to change the scopes of the files (examples can be found in existing subdirectories) and modify the CMakeLists in <strong>/libs</strong> to include the newly created subdirectory. Project scope classes should be added in the <strong>/libs</strong> folder.</p>
<p>Important definitions (variables types, constants, precision) are in the file called <strong>stemtypes_fftw3.hpp</strong>, also in the <strong>/libs</strong> folder. New definitions should be appended to this file.</p>
<p><strong>Parallelization with CUDA/ArrayFire</strong></p>
<p>ArrayFire will be imported automatically when building the project with cmake.</p>
<p>When creating CUDA source files (files with .cu extensions), make sure it stands alone and does not include any non-native C++ libraries. The NVCC compiler is known to have compilation issues with fftw3, Boost, etc. Additionally, modify the CMakeList.txt in the same directory to explicitly compile .cu files. An example of implementing raw CUDA files is <strong>CUDA2DPotential.cu</strong> in <strong>/libs/potentials</strong>. Using ArrayFire along side CUDA memory allocation (cudaMalloc()) and libraries (cuBLAS, cuFFT, etc.) sometimes produce wrong results or even memory violations. Therefore, implement CUDA kernels only when absolutely needed as most operations have equivalent ArrayFire counterparts.</p>
<p><strong>Main Executable</strong></p>
<p>The main <strong>stem3.cpp</strong> is in the <strong>/stem3</strong> folder. The <strong>Bootstrapper</strong> class handles creation of wavefunctions, potentials, etc. Append to the appropriate <strong>Register[TYPENAME]Types()</strong> functions when implementing new features.</p>
<p><strong>New Config Parameters</strong></p>
<p>To add new external parameters, modify appropriate sections of the <strong>read_qsc</strong> class in <strong>/libs/Config_IO</strong>.</p>
</div>
<div class="section" id="license">
<h1>License<a class="headerlink" href="#license" title="Permalink to this headline">¶</a></h1>
<p>Slicepp is licensed under the GNU General Public License v2.0, please see License.txt for more details.</p>
</div>

  <h3><a href="#">Table Of Contents</a></h3>
  <ul>
<li><a class="reference internal" href="#about">About</a></li>
<li><a class="reference internal" href="#installation">Installation</a></li>
<li><a class="reference internal" href="#execution">Execution</a></li>
<li><a class="reference internal" href="#programming-guide">Programming Guide</a></li>
<li><a class="reference internal" href="#license">License</a></li>
</ul>
</ul>
</div>

    

    
  </body>
</html>
