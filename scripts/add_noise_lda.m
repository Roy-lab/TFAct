% edgefile : tab delimited; two columns. Col1: regulator; col2: target
% k: lda topic number
% multiplicity_word_sampled: the number of words sampled with replacement for each regulator (document)
%   equals multiplicity_word_sampled * the number of outgoing edges from the regulator. I used 2 for this.
%   also, if an edge has been sampled twice or more, it has confidence of 1. Else 0.8.
% prefix: save lda model as ${outdir}/lda_model_k${k}_${prefix}.mat and 3-column tab separated edge list
%   as ${outdir}/${prefix}_prior.txt
% outdir: output directory

function [mdl10]=add_noise_lda(edgefile, k, multiplicity_word_sampled, prefix, outdir)

    edges = readtable(edgefile, 'FileType', 'text', 'Delimiter', '\t', 'ReadVariableNames', false);

    parents = unique(edges.Var1);
    children = unique(edges.Var2);
    m = length(parents);
    n = length(children);
    parentMap = containers.Map(parents, 1:m);
    childMap = containers.Map(children, 1:n);

    d = zeros(m, n);

    for i = 1:size(edges, 1)
        parent = string(edges.Var1(i));
        child = string(edges.Var2(i));
        r = parentMap(parent);
        c = childMap(child);
        d(r, c) = 1;
    end

    mdl10=fitlda((ceil(d*100)),k);

    save(sprintf('%s/lda_model_k%d_%s.mat',outdir,k,prefix),'mdl10','-v7.3');
    % load(sprintf('%s/lda_model_k%d_%s%s%s.mat',outdir,k,dat,prefix,cf))   

    document_topic_distr=mdl10.DocumentTopicProbabilities;
    topic_word_distr=mdl10.TopicWordProbabilities;
    num_docs = size(document_topic_distr, 1);
    num_words = size(topic_word_distr, 1);
    num_topics = size(topic_word_distr, 2);

    sim_mat = zeros(num_docs, num_words);


    for doc = 1:num_docs
        new_words=zeros(1, num_words);
        for cur_word = 1:ceil(multiplicity_word_sampled*sum(d(doc,:)))
            topic = mnrnd(1, document_topic_distr(doc, :));
            topic_index = find(topic);
            
            word_probs = topic_word_distr(:,topic_index);
            new_words = new_words +  mnrnd(1, word_probs);
        end
        
        sim_mat(doc, :) = new_words;
    end
    sim_mat2=sim_mat;
    sim_mat2(sim_mat2>0)=1;
    fprintf("for %s, num edges in original net: %d; num edges in simulated net: %d; num edges in common: %d\n", prefix, sum(sum(d)),sum(sum(sim_mat2)), sum(sum(sim_mat2.*d)));

    sim_mat(sim_mat==1)=0.8;
    sim_mat(sim_mat>1)=1;
    [row, col] = find(sim_mat);
    conf = sim_mat(sim_mat ~= 0);
    edgeList = [row, col, conf];

    namedEdgeList = cell(size(edgeList, 1), 3);
    for i = 1:size(edgeList, 1)
        namedEdgeList{i, 1} = parents{edgeList(i, 1)};
        namedEdgeList{i, 2} = children{edgeList(i, 2)};
        namedEdgeList{i, 3} = edgeList(i, 3);
    end

    outfile=sprintf("%s/%s_prior.txt", outdir, prefix);
    f = fopen(outfile, 'w');

    for i = 1:size(namedEdgeList, 1)
        fprintf(f, '%s\t%s\t%.2f\n', namedEdgeList{i, 1}, namedEdgeList{i, 2}, namedEdgeList{i, 3});
    end

    fclose(f);