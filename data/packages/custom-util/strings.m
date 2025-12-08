function result = titlecase(str)
  % Convert to lowercase first
  str = lower(str);

  % Make the first character uppercase
  if length(str) > 0
    str(1) = upper(str(1));
  endif

  % Find spaces and capitalize the letter after each space
  spaces = find(str == ' ');
  for i = 1:length(spaces)
    if spaces(i) < length(str)
      str(spaces(i) + 1) = upper(str(spaces(i) + 1));
    endif
  endfor

  result = str;
endfunction

printf('String functions loaded ✓\n');
