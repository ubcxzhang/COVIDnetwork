import java.io.*;
import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.List;
import java.util.function.Supplier;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import org.apache.commons.math3.distribution.NormalDistribution;

import it.unimi.dsi.webgraph.labelling.ArcLabelledImmutableGraph;
import it.unimi.dsi.webgraph.labelling.Label;

public class K_BZ_new_biominProject_COVID {

    static NormalDistribution norm = new NormalDistribution(); //same as before

    ArcLabelledImmutableGraph G; //same as before
    int n; //same as before
    double eta; //same as before
    int precision; //same as before
    String DPtype; //same as before
    int L; //same as before
    int m; //same as before

    double threshold_sumP;
    int threhsold_etaDegree;
    String whiteListFilename; // filename for list of node to skip screening

    int processors; //same as before

    int[] etadeg; //same as before
    BitSet gone; //same as before
    BitSet invalid; //same as before
    BitSet whiteList; // list of node to skip in screening

    int[] deg; // array of lower bounds //same as before
    int[] pos; //same as before


    //step 1: deleting vertices based on their value of sum(pi)
    //for a vertex v, deleted_step1(v) returns true if either sum(pi) < threshold_sumP
    //and deleted_step1(v) returns false if sum(pi) >= threshold_sumP
    //in fact, if deleted_step1(v) is true, it means that v has been deleted.
    final BitSet deleted_step1;


    //for those vertices that survive from step 1, deleted_step2(v) returns true if \eta-degree of v (obtined by CLT) is less than threhsold_etaDegree
    //and deleted_step2(v) returns false if \eta-degree of v is greater than or equal to threshold.
    //in fact, if deleted_step2(v) is true, it means that v has been deleted.
    final BitSet deleted_step2;

    int num_alive; //stores number of vertices


    int[] bin; //same as before

    int[] alive_vertices; //stores the index of vertices that are alive after step 1 and step 2
    int[] vert; //same as before

    int md; // max eta-degree

    public K_BZ_new_biominProject_COVID(String basename, double eta, int L, int precision, double threshold_sumP,
                                        int threhsold_etaDegree, String whiteListFilename, String DPtype) throws Exception {

        G = ArcLabelledImmutableGraph.load(basename);
        m = (int) (G.numArcs() / 2);
        n = G.numNodes();

        this.eta = eta;
        this.DPtype = DPtype;
        this.L = L;
        this.precision = precision;

        this.threshold_sumP = threshold_sumP;
        this.threhsold_etaDegree = threhsold_etaDegree;
        this.whiteListFilename = whiteListFilename;

        // processors = Runtime.getRuntime().availableProcessors();
        // align num of cpu with compute canada job script
        processors = 7;
        System.out.println("gy: Number of processors is " + processors);

        etadeg = new int[n];
        gone = new BitSet(n);
        invalid = new BitSet(n);
        whiteList = new BitSet(n); // a bitset with length = full node number in graph

        readWhiteList(whiteListFilename, whiteList);

        pos = new int[n];
        deg = new int[n];
        Arrays.fill(deg, -1); //new: -1 is for vertices removed
        Arrays.fill(pos, -1); //new: -1 is for vertices removed

        // it keeps track of vertices which are removed based on their sumP values
        deleted_step1 = new BitSet(n);

        long startTime = System.currentTimeMillis();


        remove_byFilter_sumP_Sequential(); //the same function, but run sequentially

        deleted_step2 = (BitSet) deleted_step1.clone();
        //just for test
        //System.out.println(deleted_step1.cardinality() + " " + deleted_step2.cardinality());

//        this.num_alive = compute_etaDegree_foraliveVertices_filter_Parallel();
        //System.out.println(n-deleted_step2.cardinality());
        this.num_alive = compute_etaDegree_foraliveVertices_filter_Sequential();

        System.out.println("gy: Total number of alive vertices are:" + num_alive + " check if this is the same with cnt_alive?");

        //now store them in an array
        alive_vertices = IntStream.range(0, deg.length).filter(i -> (deg[i] >=threhsold_etaDegree)).toArray(); // for all node i that has deg[i]>= threhsold_etaDegree, store i in an array
        System.out.println("gy: alive_vertices array length (should be cnt_alive):" + alive_vertices.length);
        //now compute lower-bound for alive vertices (same as before in peeling process)

        md = 0;



        for (int i=0;i<num_alive;i++) {
            int u=alive_vertices[i];
            if (deleted_step2.get(u) == false) { //the vertex has not been removed in the preprocess stage

                int temp = LYc(G.labelArray(u), G.successorArray(u), deleted_step2) - (eta > 0 ? 1 : 0);
                // We correct by -1 because, in very few cases, LYc overestimates and doesn't
                // produce a lower bound.
                // We only do the correction if eta>0.
                // If eta=0 (treating the graph as deterministic), LYc is very accurate, so no need to correct.

                deg[u] = temp < 0 ? 0 : temp;




                if (deg[u] >= md) {
                    md = deg[u];
                }

                invalid.set(u, true);
            }
        }

        bin = new int[md + 1]; // md+1 because we can have zero degree
        vert = new int[num_alive]; //keep edges sorted based on their eta-degree

        long endTime = System.currentTimeMillis();

        System.out.println("Time elapsed for initial (sec) = " +(endTime - startTime) / 1000.0);

        System.out.println("n=" + n+  "," + "m=" + m + "," + " Max_eta-degree=" + md);



        /*
         * IntStream.range(0, processors) .parallel() .forEach(processor -> {
         * computeInitialEtaDegs_inProcessor(processor); });
         */

    }

    public void writeTest(int[] arr, String filename) throws IOException {
        BufferedWriter w = new BufferedWriter(new FileWriter(filename));

        List<Integer> l = new ArrayList<>();
        for (int i = 0; i < arr.length; i++) {
            int u = arr[i];
            l.add(u);

        }

        for (int i=0;i<l.size();i++) {
            int u = l.get(i);
            w.write(u);
            w.write("\n");
        }
        w.close();
    }

    public int compute_etaDegree_foraliveVertices_filter_Parallel() {

        BigInteger cnt_alive = IntStream.range(0, n).parallel().map(u -> {

            if (deleted_step1.get(u) == false) { //if u is still alive
                ArcLabelledImmutableGraph H = G.copy();
                BitSet deleted_copy = (BitSet) deleted_step1.clone();
                Label[] u_label = H.labelArray(u);
                int[] neighborArray = H.successorArray(u);

                int temp = LYc(u_label, neighborArray, deleted_copy);
                deg[u] = temp < 0 ? 0 : temp;



                if (deg[u] >= threhsold_etaDegree) {
                    return 1;
                } else { //the vertex should be removed
                    deleted_step2.set(u);
                    return 0;
                }
            }else {
                return 0;
            }
        }).mapToObj(BigInteger::valueOf).reduce(BigInteger.ZERO, BigInteger::add);

        System.out.println("Number of vertices alive after filtering by eta-degree (cnt_alive)=" + cnt_alive + "==" + (n-deleted_step2.cardinality()));

        return cnt_alive.intValue();


    }

    public int compute_etaDegree_foraliveVertices_filter_Sequential() {
        // gy: already modified for COVID analysis, cannot smoothly change to parallel ver above!


        int cnt_alive=0; // alived after CLT + nodes in white list

        for (int u = 0; u < n; u++) {
            // if current node u is in white list, skip calculation
            if (whiteList.get(u)) {
                deg[u] = threhsold_etaDegree; // set the skipped node's CLT result to be the threshold to prevent it from deletion
                cnt_alive++;
                continue;
            }
            if (deleted_step1.get(u) == false) { //if u is still alive
                Label[] u_label = G.labelArray(u);
                int[] neighborArray = G.successorArray(u);


                int temp = LYc(u_label, neighborArray, deleted_step1);
                deg[u] = temp < 0 ? 0 : temp;
                if (deg[u] >= threhsold_etaDegree) {
                    cnt_alive++;
                }else {
                    deleted_step2.set(u); // u is deleted in step 2
                }

            }

        }

        System.out.println("Number of vertices alive after filtering by eta-degree (cnt_alive)=" + cnt_alive + "==" + (n-deleted_step2.cardinality()));
        return cnt_alive;

    }

    public int LYc(Label[] LabelArray, int[] neighborArray, BitSet deleted_copy) {

        double s_n2 = 0.0;
        double expectedValue = 0.0;
        double s_n;
        double z;
        int core = 0;

        int count_neigh_alive = 0;

        // the for-loop computes sn^2 and the sum of mu_i
        for (int j = 0; j < LabelArray.length; j++) {

            int v = neighborArray[j]; // neighbor vertex

            // make sure correct for parallel
            if (deleted_copy.get(v) == false) { // if v has not been removed in the previous stage
                double prob_value = LabelArray[j].getLong() * (1.0 / Math.pow(10, precision));
                s_n2 = s_n2 + (prob_value * (1 - prob_value));
                expectedValue = expectedValue + (prob_value);
                count_neigh_alive++;
            }

        }

        if (count_neigh_alive == 0) {
            return 0;
        } else {
            s_n = Math.sqrt(s_n2);

            double prob;

            for (int k = 1; k <= count_neigh_alive; k++) {
                z = (double) ((k - 0.5) - expectedValue) / s_n;
                prob = 1 - norm.cumulativeProbability(z);

                if (prob < eta) {
                    core = k - 1;
                    break;
                }
                if (k == count_neigh_alive) {
                    core = count_neigh_alive;
                }

            }
            return core;
        }
    }

    public void remove_byFilter_sumP_Sequential() {


        int count = 0;
        for (int u = 0; u < n; u++) {
            Label[] u_label = G.labelArray(u);
            double[] sumP_sum = find_sumP_sum(u_label);
            // if current node u is in white list, skip calculation
            if (whiteList.get(u)) {
                count++;
                continue;
            }
            if ((sumP_sum[0] < threshold_sumP)) {
                deleted_step1.set(u);

            } else {
                count++;
            }
        }
        System.out.println(
                "Number of vertices alive after filtering by sumP_sum=" + count + "==" + (n - deleted_step1.cardinality()));
        System.out.println("Vertices_deleted=" + deleted_step1.cardinality() + " n=" + n);
    }

    public double[] find_sumP_sum(Label[] u_label) {
        double[] sumP_sum = new double[2];

        double sum_for_p = 0.0;

        //prob is used to change the long value of an edge to its actual probability.
        //for instance, when precision=2 in a webgraph, for an edge e with probability 0.8,
        //we store this value as 80. Then, prob is equal to 0.01. By multiplying 80*(0.01)
        //we get the probability of the edge
        double prob = 1.0 / Math.pow(10, precision); //p_i=0.8 in webgraph we store it as 80


        for (int j = 0; j < u_label.length; j++) {

            Long val = u_label[j].getLong();

            double sumP = val * prob; //sumP > temP

            sum_for_p += sumP;

        }
        sumP_sum[0] = sum_for_p; //sum(pi)

        return sumP_sum;
    }

    private void computeInitialEtaDegs_inProcessor(int processor) {
        // computing range of nodes for which to compute initial eta-deg
        int[] range = getRangeOfNodes(processor);

        ArcLabelledImmutableGraph H = G.copy();

        for (int v = range[0]; v < range[1]; v++) {

            int temp = KCoins.LYc(H.labelArray(v), H.outdegree(v), H.outdegree(v), eta, precision) - (eta > 0 ? 1 : 0);
            // We correct by -1 because, in very few cases, LYc overestimates and doesn't
            // produce a lower bound.
            // We only do the correction if eta>0.
            // If eta=0, LYc is very accurate, so no need to correct.

            deg[v] = temp < 0 ? 0 : temp;

            /*
             * deg[v] = KCoins.LYc( H.labelArray(v), H.outdegree(v), H.outdegree(v),
             * eta,precision);
             */

        }
    }

    // Return an array "range" of size 2. range[0] is "from", range[1] is "to".
    private int[] getRangeOfNodes(int processor) {
        int num_nodes = n / processors;
        int[] range = new int[2];
        range[0] = processor * num_nodes;
        range[1] = processor != processors - 1 ? (processor + 1) * num_nodes : n;
        return range;
    }

    public int[] KCoreCompute() {

        for (int d = 0; d <= md; d++)
            bin[d] = 0;
        for (int i = 0; i < alive_vertices.length; i++) { //new
            int v = alive_vertices[i]; //vertex at position i
            bin[deg[v]]++;
        }
		
		/*int count=0;
		for (int i=0;i<deg.length;i++) {
			if (deg[i] >= 0 && deleted_step2.get(i)==false) {
				count++;
			}
		}
		System.out.println(count + " " + count);*/

        //System.out.println("num_alive " + num_alive + " " + alive_vertices.length);

        int start = 0; // start=1 in original, but no problem (same as
        // deterministic case)
        for (int d = 0; d <= md; d++) {
            int num = bin[d];
            bin[d] = start;
            start += num;
        }

        // bin-sort vertices by degree (same as deterministic case)
        for (int i = 0; i < alive_vertices.length; i++) { //new
            int v = alive_vertices[i];
            pos[v] = bin[deg[v]];
            vert[pos[v]] = v;
            bin[deg[v]]++;
        }
        // recover bin[] (same as deterministic case)
        for (int d = md; d >= 1; d--)
            bin[d] = bin[d - 1];
        bin[0] = 0; // 1 in original (same as deterministic case)


        //for test
        //for (int i=0;i<alive_vertices.length;i++) {
        //System.out.println(vert[i] + " " + deg[vert[i]]);
        //}



        // main algorithm
        for (int i = 0; i < num_alive;) {

            int v = vert[i]; // smallest degree vertex that will be deleted

            if (invalid.get(v) == false) {

                gone.set(v);

                int v_deg = G.outdegree(v);

                int[] N_v = G.successorArray(v);
                for (int j = 0; j < v_deg; j++) {
                    int u = N_v[j];

                    if (deleted_step2.get(u) == false) { //if u has not been deleted

                        if (deg[u] == deg[v] && invalid.get(u) == true) {

                            recompute_and_swap(u);
                            invalid.clear(u);
                        }

                        if (deg[u] > deg[v]) {
                            swap_left(u);
                            invalid.set(u);
                        }

                    }
                }
                i++;
            } else {
                recompute_and_swap(v);
            }
        }


        return deg;
    }

    void recompute_and_swap(int v) {
        etadeg[v] = recompute(v);
        invalid.clear(v);
        int diff = etadeg[v] - deg[v];

        for (int d = 0; d < diff; d++)
            swap_right(v);
    }

    void swap_left(int u) {
        int du = deg[u];
        int pu = pos[u];
        int pw = bin[du];
        int w = vert[pw];
        if (u != w) {
            pos[u] = pw;
            vert[pu] = w;
            pos[w] = pu;
            vert[pw] = u;
        }
        bin[du]++;
        deg[u]--;

    }

    void swap_right(int u) {
        int du = deg[u];
        int pu = pos[u];
        int pw = bin[du + 1] - 1;
        int w = vert[pw];
        if (u != w) {
            pos[u] = pw;
            vert[pu] = w;
            pos[w] = pu;
            vert[pw] = u;
        }
        bin[du + 1]--;
        deg[u]++;
    }

    int recompute(int v) {
        int v_deg_real;

        // rebuild label_array for v
        v_deg_real = 0; // real degree, i.e. num of neighbors still alive
        int[] v_neighbors = G.successorArray(v);
        int v_deg = G.outdegree(v);
        for (int t = 0; t < v_deg; t++) {
            int u = v_neighbors[t];

            if ((gone.get(u) == false) && (deleted_step2.get(u) == false))
                v_deg_real++;
        }

        int Index = 0;

        Label[] v_label_real = new Label[v_deg_real];
        Label[] v_label = G.labelArray(v);
        for (int t = 0; t < v_deg; t++) {
            int u = v_neighbors[t];

            if (gone.get(u) == false && (deleted_step2.get(u) == false)) {
                v_label_real[Index] = v_label[t];
                Index++;
            }
        }

        if (v_deg_real > L) {
            return LYc(v_label_real, v_neighbors, deleted_step2);
        } else {
            //if (DPtype.equals("DP_log2"))
            //return KCoins.DP_log2(v_label_real, v_deg_real, v_deg_real, eta, precision);
            //else
            return DP(v_label_real, v_deg_real, v_deg_real);
        }

    }

    public int DP(Label[] u_labels, int H, int J) {

        double etaProbability = 1.0;
        double[][] dp = new double[H + 1][2];
        int sw = 0;

        for (int j = 0; j <= J; j++) { // column

            for (int h = 0; h <= H; h++) { // row
                if (j == 0 && h == 0) {
                    dp[h][sw] = 1;
                } else if (h < j) {
                    dp[h][sw] = 0;
                } else {

                    double prob = (1.0 / Math.pow(10, precision)) * u_labels[h - 1].getLong();
                    dp[h][sw] = prob * (j == 0 ? 0 : dp[h - 1][1 - sw]) + (1- prob) * dp[h - 1][sw];
                }
            }

            double probability = dp[H][sw];


            etaProbability -= probability;

            if (etaProbability < eta) {
                return j;
            }
            sw = 1 - sw;
        }

        if (J == H) {
            return H;
        } else {
            return J;
        }

    }

    // read in whiteList File and set the corresponding cell in whiteList according to nodes in File
    public void readWhiteList(String filename, BitSet whiteList) throws IOException {
        int counter = 0;
        int u;

        BufferedReader br = new BufferedReader(new FileReader(filename));
        String line = br.readLine(); // whiteList file will be one column with each line being a `integer node name`
        while (line != null) {
//            u = Integer.parseInt(line.split(";")[0]); // if each line is `integer node name;`
            u = Integer.parseInt(line);
            whiteList.set(u);
            counter++;
            line = br.readLine();
        }
        System.out.println("gy: total number of nodes to skip screening:" + counter);
    }

    public void writeResults(int[] core, String filename) throws IOException {
        BufferedWriter w = new BufferedWriter(new FileWriter(filename));
        for (int v = 0; v < n; v++) {
            if (deleted_step2.get(v) == false) {
                w.write(v + "\t" + core[v]);
                w.write("\n");
            }
        }
        w.close();
    }

    public static void main(String[] args) throws Exception {

        //args = new String[] { "biomine-proc.w", "0.1", "1500", "16", "20", "20" };
        //gy 210722 omitted for BIB paper added experiments
//        args = new String[] { "Biomine_data_unique_withCOVID_renamed-proc.w", "0.5", "1000", "17", "5", "10", "node_to_skip.txt"};


        if (args.length < 7) {
            System.err.println("Specify: basename eta L precision thresh_sumPQ threhsold_etaDegree WhiteList_filename DPtype(optional)");
            System.exit(1);
        }

        System.out.println("Basename: " + args[0] + " eta: " + args[1] + " L: " + args[2] +" precision: " + args[3] + " stage-1-threshold: " + args[4] + " stage-2-threshold: " + args[5] + "WhiteList_filename: " + args[6]);


        String basename = args[0];
        double eta = Double.parseDouble(args[1]);
        int L = Integer.parseInt(args[2]);
        int precision = Integer.parseInt(args[3]);
        double thresh_sumPQ = Double.parseDouble(args[4]);
        int threhsold_etaDegree = Integer.parseInt(args[5]);
        String whiteListFilename = args[6];
        String DPtype = args.length == 8 ? args[7] : "DP";

        System.out.println("Starting " + basename);
        long startTime = System.currentTimeMillis();

        K_BZ_new_biominProject_COVID tt = new K_BZ_new_biominProject_COVID(basename, eta, L, precision, thresh_sumPQ,
                threhsold_etaDegree, whiteListFilename, DPtype);

        //new line: save all the nodes (based on their index) that survive

        int num_nodes = tt.n;
        String filename = "AliveNodes"+basename+"-"+eta+".txt";
        BufferedWriter w = new BufferedWriter(new FileWriter(filename));
        w.write("Index of nodes that are alive");
        w.write("\n");

        for (int i=0;i<num_nodes;i++) {
            if (tt.deleted_step2.get(i)==false) { //if node has not been deleted at the end of filtering
                w.write(i + "\n");
            }
        }
        w.close();


        int[] res = tt.KCoreCompute();
        System.out.println(args[0] + ": Time elapsed (sec) = " + (System.currentTimeMillis() - startTime) /1000.0);

        int kmax = 0;
        double sum = 0;
        int cnt = 0;
        for (int i=0;i<res.length;i++) {
            if (tt.deleted_step2.get(i) == false) {
                //if (res[i] >=0) {
                if(res[i] >=kmax) {
                    kmax=res[i];
                }
                sum += res[i];
                cnt++;
                //}
            }

            //if (res[i] > 0) {
            //	 cnt++;
            //}
        }
        System.out.println("|V|	|E|	etadmax	kmax	kavg");
        System.out.println(cnt + "\t" + tt.m + "\t" + tt.md + "\t" + kmax + "\t" + (sum / cnt));
        System.out.println();

        //it stores the core number of the alive vertices only
        tt.writeResults(res, basename+"eta-" + eta + "-bz.txt");






    }
}