function txt = center_text(txt, width, within, add_space)
if nargin<3
    within = ' ';
end
if nargin<4
    add_space = false;
end

if add_space
    txt = [' ' txt ' '];
end

l = length(txt);

m = width - l;
m_left = floor(m/2);
m_right = m - m_left;

txt = [repmat(within, 1, m_left) txt repmat(within, 1, m_right)];

end
