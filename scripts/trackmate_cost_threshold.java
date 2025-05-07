// Custom TrackMate linking-cost threshold snippet
// mCostThreshold = MaxDistance^2
final double cost = mCostFunction.linkingCost(source, target);
// Retrieve the "RADIUS" feature from source and target spots
final double sradius = source.getFeature("RADIUS");
final double tradius = target.getFeature("RADIUS");
// If sqrt(cost) exceeds combined radii + sqrt(mCostThreshold), skip linking
if (Math.sqrt(cost) > (sradius + tradius + Math.sqrt(mCostThreshold))) {
    continue;  // Do not link these spots
}
