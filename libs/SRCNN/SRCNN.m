function im_h = SRCNN(model, im_b)

%% load CNN model parameters
load(model);
[conv1_patchsize2,conv1_filters] = size(weights_conv1);
conv1_patchsize = sqrt(conv1_patchsize2);
[conv2_filters,conv3_patchsize2] = size(weights_conv3);
conv3_patchsize = sqrt(conv3_patchsize2);
[hei, wid] = size(im_b);

%% conv1
weights_conv1 = reshape(weights_conv1, conv1_patchsize, conv1_patchsize, conv1_filters);
conv1_data = zeros(hei, wid, conv1_filters);
for i = 1 : conv1_filters
    conv1_data(:,:,i) = imfilter(im_b, weights_conv1(:,:,i), 'same', 'replicate');
    conv1_data(:,:,i) = max(conv1_data(:,:,i) + biases_conv1(i), 0);
end

%% conv2
conv2_data = zeros(hei, wid, conv2_filters);
for i = 1 : conv2_filters
    for j = 1 : conv1_filters
        conv2_data(:,:,i) = conv2_data(:,:,i) + weights_conv2(j,:,i) * conv1_data(:,:,j);
    end
    conv2_data(:,:,i) = max(conv2_data(:,:,i) + biases_conv2(i), 0);
end

%% conv3
conv3_data = zeros(hei, wid);
for i = 1 : conv2_filters
    conv3_subfilter = reshape(weights_conv3(i,:), conv3_patchsize, conv3_patchsize);
    conv3_data(:,:) = conv3_data(:,:) + imfilter(conv2_data(:,:,i), conv3_subfilter, 'same', 'replicate');
end

%% SRCNN reconstruction
im_h = conv3_data(:,:) + biases_conv3;