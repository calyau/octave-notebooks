function ensure_dirs(filepath)
  % Extract directory from filepath
  [dir, name, ext] = fileparts(filepath);

  % Create directory if it doesn't exist
  if !exist(dir, 'dir')
    mkdir(dir)
  endif
endfunction

printf('File system functions loaded ✓\n');
