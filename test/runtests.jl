using Test, Brunnermeier_Sannikov_Replication

@test hello("Julia") == "Hello, Julia"
@test domath(2.0) ≈ 7.0
