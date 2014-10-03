import java.io.BufferedReader;
import java.io.FileReader;
import java.util.Arrays;


public class RunTests {
	public static String[] TEST_FILES = {
		"a280_n279_bounded-strongly-corr_01",
		"a280_n1395_uncorr-similar-weights_05",
		"a280_n2790_uncorr_10",
		"fnl4461_n4460_bounded-strongly-corr_01",
		"fnl4461_n22300_uncorr-similar-weights_05",
		"fnl4461_n44600_uncorr_10",
		"pla33810_n33809_bounded-strongly-corr_01",
		"pla33810_n169045_uncorr-similar-weights_05",
		"pla33810_n338090_uncorr_10"
	};
	public static String[] ALGORITHMS = {
		"alg1",
		"alg2"
	};
	public static void main(String[] args) {
		int numIterations = 3;
		int numTimeSaves = 3;
		// Run all tests
		for (int algNum = 0; algNum < ALGORITHMS.length; algNum++) {
			for (int testNum = 0; testNum < TEST_FILES.length; testNum++) {
				for (int iteration = 0; iteration < numIterations; iteration++) {
					Driver.runTest("instances/"+TEST_FILES[testNum], algNum, iteration);
				}
			}
		}
		
		// Combine results
		
		for (int algNum = 0; algNum < ALGORITHMS.length; algNum++) {
			for (int testNum = 0; testNum < TEST_FILES.length; testNum++) {
				double[][] obj = new double[numIterations][numTimeSaves];
				for (int iteration = 0; iteration < numIterations; iteration++) {
					String resultFile = String.format("%s_%s_%d",
							ALGORITHMS[algNum],
							TEST_FILES[testNum],
							iteration);
					BufferedReader br = new BufferedReader(new FileReader(resultFile));
					for (int timeIdx = 0; timeIdx < obj[iteration].length; timeIdx++) {
						obj[iteration][timeIdx] = Integer.parseInt(br.readLine());
					}
				}
				// Calc average, standev, best (include iteration num that was best)
				double[] objAvg = new double[numTimeSaves];
				double[] objBest = new double[numTimeSaves];
				double[] objBestIter = new double[numTimeSaves];
				Arrays.fill(objAvg, 0);
				Arrays.fill(objBest, Double.NEGATIVE_INFINITY);
				Arrays.fill(objBestIter, -1);
				for (int iteration = 0; iteration < numIterations; iteration++) {
					for (int timeIdx = 0; timeIdx < objAvg.length; timeIdx++) {
						objAvg[timeIdx] += obj[iteration][timeIdx];
						if (obj[iteration][timeIdx] > objBest[timeIdx]) {
							objBest[timeIdx] = obj[iteration][timeIdx];
						}
					}
				}
			}
		}
	}
}
