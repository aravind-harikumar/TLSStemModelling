function visulizelayerdata(I)

for i = 1:64
    subplot(8,8,i)
    imagesc(I(:,:,i))
end

end